# ---------------------------------------------------------------------------
# check_figure() — verify an EPS/PS figure against Taylor & Francis / JASA specs
#
# Checks, all performed on the EPS ITSELF (never on a ghostscript conversion,
# because gs re-embeds fonts and can turn vector glyphs into Type 3 bitmaps —
# both of which corrupt the verdict):
#   1. file really is EPSF
#   2. page size from %%BoundingBox, vs intended size
#   3. every font is SUPPLIED (embedded), none merely NEEDED
#   4. no font outside the publisher's standard list
#   5. no Type 3 (bitmap-glyph) fonts
#   6. no raster images (figure is pure vector)
#   7. minimum stroke width, measured in real points, after any reduction
#
# Optional cross-check with ghostscript: if gs must load a font of its own to
# interpret the file, that font was not embedded.
#
# Zero package dependencies (base R only).
# ---------------------------------------------------------------------------

STANDARD_FONTS <- c("Times", "Times-Roman", "Times-Bold", "Times-Italic",
                    "Times-BoldItalic", "TimesNewRoman", "TimesNewRomanPSMT",
                    "Helvetica", "Helvetica-Bold", "Helvetica-Oblique",
                    "Helvetica-BoldOblique", "Arial", "Arial-Bold",
                    "ArialMT", "Arial-BoldMT", "Symbol",
                    "Courier", "Courier-Bold", "Courier-Oblique",
                    "Courier-BoldOblique")

check_figure <- function(path,
                         width_in    = NULL,   # intended device width  (inches)
                         height_in   = NULL,   # intended device height (inches)
                         min_line_pt = 0.5,    # 0.5 = T&F/CRC; 0.3 = T&F journals
                         scale       = 1,      # publisher reduction (0.85 = 15% down)
                         min_text_pt = 8,      # ASA minimum at final printed size
                         base_pt     = NULL,   # the pointsize/base_size you used
                         use_gs      = TRUE) {

  stopifnot(file.exists(path))
  raw <- readBin(path, "raw", file.size(path))
  txt <- rawToChar(raw); Encoding(txt) <- "latin1"
  lines <- strsplit(txt, "\r\n|\n|\r")[[1]]

  res <- list(); ok <- TRUE
  add <- function(name, pass, detail) {
    res[[length(res) + 1]] <<- list(name = name, pass = pass, detail = detail)
    if (!isTRUE(pass)) ok <<- FALSE
  }

  ## 1. EPSF conformance ----------------------------------------------------
  hdr <- lines[1]
  is_eps <- grepl("^%!PS-Adobe.*EPSF", hdr)
  add("EPSF header", is_eps, hdr)

  ## 2. page size -----------------------------------------------------------
  bb <- grep("^%%BoundingBox:", lines, value = TRUE)
  bb <- bb[!grepl("atend", bb, ignore.case = TRUE)][1]
  if (is.na(bb)) {
    add("BoundingBox", FALSE, "absent")
    w_in <- h_in <- NA
  } else {
    v <- as.numeric(strsplit(sub("^%%BoundingBox:\\s*", "", bb), "\\s+")[[1]])
    w_pt <- v[3] - v[1]; h_pt <- v[4] - v[2]
    w_in <- w_pt / 72;   h_in <- h_pt / 72
    det <- sprintf("%.0f x %.0f pt = %.3f x %.3f in", w_pt, h_pt, w_in, h_in)
    pass <- TRUE
    if (!is.null(width_in))
      pass <- pass && abs(w_in - width_in) < 0.02
    if (!is.null(height_in))
      pass <- pass && abs(h_in - height_in) < 0.02
    if (!is.null(width_in))
      det <- paste0(det, sprintf("  (intended %.2f x %.2f in)", width_in,
                                 if (is.null(height_in)) NA else height_in))
    add("Page size", pass, det)
  }

  ## 3/4/5. fonts -----------------------------------------------------------
  # NEEDED = the interpreter must supply it = NOT embedded.
  needed <- character()
  ndx <- grep("^%%DocumentNeededResources:", lines)
  if (length(ndx)) {
    i <- ndx[1]
    needed <- c(needed, sub(".*font\\s+", "", lines[i]))
    j <- i + 1
    while (j <= length(lines) && grepl("^%%\\+", lines[j])) {
      if (grepl("font", lines[j])) needed <- c(needed, sub(".*font\\s+", "", lines[j]))
      j <- j + 1
    }
  }
  supplied <- sub("^%%BeginResource:\\s*font\\s+", "",
                  grep("^%%BeginResource:\\s*font", lines, value = TRUE))
  needed <- unique(trimws(needed)); supplied <- unique(trimws(supplied))

  add("Fonts embedded",
      length(needed) == 0,
      if (length(needed))
        paste0("NOT embedded: ", paste(needed, collapse = ", "),
               "  |  embedded: ", if (length(supplied)) paste(supplied, collapse = ", ") else "none")
      else paste0("all embedded: ",
                  if (length(supplied)) paste(supplied, collapse = ", ") else "(no text)"))

  allf <- unique(c(needed, supplied))
  base_name <- sub("^[A-Z]{6}\\+", "", allf)          # strip subset tag
  nonstd <- base_name[!base_name %in% STANDARD_FONTS]
  add("Standard fonts only", length(nonstd) == 0,
      if (length(nonstd)) paste("non-standard:", paste(nonstd, collapse = ", "))
      else "ok")

  ft <- as.integer(sub(".*/FontType\\s+([0-9]+).*", "\\1",
                       grep("/FontType\\s+[0-9]", lines, value = TRUE)))
  ft <- ft[!is.na(ft)]
  add("No Type 3 (bitmap) fonts", !any(ft == 3),
      if (length(ft)) paste0("FontType(s) present: ",
                             paste(sort(unique(ft)), collapse = ", "),
                             "  (42 = TrueType-embedded, 1 = Type 1)")
      else "no font programs")

  ## 6. raster --------------------------------------------------------------
  # count INVOCATIONS only. Definition lines start with '/'. Hex font payload
  # can never match these words, so /sfnts data is harmless here.
  is_def <- grepl("^\\s*/", lines)
  inv <- grepl("\\b(cairo_image|cairo_imagemask)\\b", lines) & !is_def
  inv <- inv | (grepl("\\b(colorimage|imagemask)\\b|[^a-zA-Z]image\\b", lines) &
                  !is_def & !grepl("%", lines))
  n_raster <- sum(inv)
  add("Pure vector (no raster)", n_raster == 0,
      sprintf("%d raster invocation(s)%s", n_raster,
              if (n_raster) " -- usually caused by alpha transparency" else ""))

  ## 7. line widths ---------------------------------------------------------
  m <- regmatches(lines, regexec("^\\s*([0-9]*\\.?[0-9]+)\\s+(w|setlinewidth)\\s*$", lines))
  lw <- suppressWarnings(as.numeric(vapply(m, function(x) if (length(x) == 3) x[2] else NA_character_, "")))
  lw <- lw[!is.na(lw) & lw > 0]
  if (length(lw)) {
    printed <- min(lw) * scale
    add("Minimum line width",
        printed >= min_line_pt - 1e-9,
        sprintf("thinnest = %.3f pt in file; x scale %.2f -> %.3f pt printed (min %.2f). widths: %s",
                min(lw), scale, printed, min_line_pt,
                paste(sprintf("%.3f", sort(unique(lw))), collapse = ", ")))
  } else {
    add("Minimum line width", NA, "no stroked lines found")
  }

  ## 8. text size (arithmetic, not measured) --------------------------------
  if (!is.null(base_pt)) {
    printed <- base_pt * scale
    add("Text size after reduction", printed >= min_text_pt - 1e-9,
        sprintf("%.2f pt x %.2f = %.2f pt printed (min %g)", base_pt, scale, printed, min_text_pt))
  }

  ## 9. ghostscript cross-check --------------------------------------------
  gs <- Sys.which("gs")
  if (use_gs && nzchar(gs)) {
    out <- suppressWarnings(system2(gs,
      c("-dNOPAUSE", "-dBATCH", "-dNODISPLAY", "-sDEVICE=nullpage", shQuote(path)),
      stdout = TRUE, stderr = TRUE))
    ext <- unique(sub("^Loading ", "", regmatches(out, regexpr("Loading [A-Za-z0-9-]+", out))))
    add("gs needs no external font", length(ext) == 0,
        if (length(ext)) paste("gs had to load:", paste(ext, collapse = ", "))
        else "gs rendered with only the file's own fonts")
  }

  ## report -----------------------------------------------------------------
  cat(sprintf("\n=== check_figure: %s ===\n", basename(path)))
  for (r in res)
    cat(sprintf("  [%s] %-28s %s\n",
                if (is.na(r$pass)) "?" else if (r$pass) "PASS" else "FAIL",
                r$name, r$detail))
  cat(sprintf("  ---> %s\n", if (ok) "OK" else "NOT COMPLIANT"))
  invisible(list(ok = ok, checks = res, width_in = w_in, height_in = h_in))
}

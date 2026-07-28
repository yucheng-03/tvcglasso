# ---------------------------------------------------------------------------
# analysis/R/00_aes.R  --  THE single visual identity for every paper figure.
#
# Nothing else in analysis/ may define a colour, a line type, a device or a
# figure size. If a figure needs a new visual element, add it HERE so the two
# figures cannot drift apart (they did drift in the prototype scripts: the same
# hex meant CGLasso in one figure and tvmgm in another).
#
# ---------------------------------------------------------------------------
# Taylor & Francis / JASA artwork rules encoded below, with the source of each:
#
#  * "Any figure containing text ... supplied as EPS or PS, not raster"; and
#    "Please do not supply files in PDF format as these are 'locked' files and
#    incompatible with our workflow software."  -> EPS is the submission format.
#    (T&F, Submission of electronic artwork.)
#  * "All fonts should be embedded ... Any fonts that are not embedded will be
#    replaced by Courier."                       -> embed_eps() below.
#  * "Use standard fonts such as Times, Helvetica, Arial, and Symbol."
#    Empirically, JASA figures are SANS-SERIF even where the body text is serif
#    (checked in Peng et al. 2009, Cao-Lin-Li COAT, Zhang-Wang-Lin CARE,
#    Yang & Peng loggle, Danaher et al. JGL) -> Helvetica.
#  * "Minimum line weight 0.3 pt for black lines on a white background."
#    R's lwd unit is 1/96 in = 0.75 pt, so lwd >= 0.4 clears it; we use >= 0.5.
#  * "Figure files must not contain layers or transparent objects."
#    -> no alpha anywhere. EPS cannot carry it and R would rasterise the page.
#  * "Do not render captions or figure titles inside the figure file."
#    -> figures carry panel identifiers only; the caption lives in the .tex.
#  * "Figures print in black and white but appear in colour online ... do not
#    encode series identity by hue alone."       -> colour AND line type/fill.
#  * ASA Style Guide: lettering "no smaller than 8 points" at final printed
#    size; measured JASA in-figure text sits at ~7.4 pt (= caption size).
#  * Measured JASA geometry: text block 516 pt = 7.17 in full width (two-column
#    layout), single column 3.50 in. Draw AT final size; do not draw large and
#    let LaTeX shrink (that is what pushes label sizes below the minimum).
# ---------------------------------------------------------------------------

## --- methods ---------------------------------------------------------------
## Order is fixed: our method first, then baselines in the order they appear in
## the paper. Used for legend order and bar order.
METHODS <- c("tvcglasso", "CGLasso", "tvmgm", "JGL")

## Display names. The folder/method string is lowercase "tvcglasso"; the paper
## writes it TVCGLasso.
METHOD_LABEL <- c(
  tvcglasso = "TVCGLasso",
  CGLasso   = "CGLasso",
  tvmgm     = "tvmgm",
  JGL       = "JGL"
)

## --- palette ---------------------------------------------------------------
## Paul Tol's "high contrast" qualitative scheme (SRON/EPS/TN/09-002), the one
## published palette designed to be simultaneously colour-blind safe AND to
## survive greyscale conversion.
##
## Verified numerically (sRGB -> CIELAB, CIEDE2000, Machado 2009 and Vienot 1999
## dichromacy models):
##   method      hex      CIE L*   greyscale hex
##   tvcglasso   #000000    0.0      #000000
##   CGLasso     #004488   29.2      #454545
##   tvmgm       #BB5566   49.5      #767676
##   JGL         #DDAA33   72.5      #B2B2B2
## Minimum adjacent L* gap 20.3; minimum pairwise CIEDE2000 after desaturation
## 17.8. All three luminance metrics (WCAG, Rec.709, L*) order the four series
## identically, so the greyscale ranking is stable whichever conversion the
## printer's RIP applies.
##
## tvcglasso is BLACK on purpose: it is the focal series, black is the only
## achromatic anchor (invariant under every colour-vision deficiency) and it
## carries the most ink in greyscale. To make it a colour instead, change the
## one line below -- everything downstream follows.
METHOD_COL <- c(
  tvcglasso = "#000000",
  CGLasso   = "#004488",
  tvmgm     = "#BB5566",
  JGL       = "#DDAA33"
)

## Redundant encoding so hue is never load-bearing (T&F greyscale rule).
## Values are R lty codes; the dash patterns are chosen to stay distinguishable
## at the ~1.8 in panel width of the 4x4 grid.
METHOD_LTY <- c(
  tvcglasso = "solid",
  CGLasso   = "22",     # short dash
  tvmgm     = "12",     # dotted
  JGL       = "4212"    # dash-dot
)

## Bars cannot carry a line type, so the LUMINANCE ordering above is the
## redundant, non-hue channel: the four fills desaturate to #000000, #454545,
## #767676, #B2B2B2 -- four evenly spaced greys spanning nearly the full ink
## range, in a fixed order that the legend states. Hatching was tried and
## rejected: it is invisible on the black series, so it distinguishes only
## three of the four bars while adding visual noise to those three.

## Point symbol, used only where a discrete operating point is marked.
METHOD_PCH <- c(tvcglasso = 21, CGLasso = 22, tvmgm = 24, JGL = 23)

## --- line weights ----------------------------------------------------------
## MEASURED, not assumed: R's `lwd` unit is 1/96 inch, so the width that lands
## in the file is lwd * 0.75 pt (verified lwd 0.3/0.4/0.5/1/2 -> 0.225/0.300/
## 0.375/0.750/1.500 pt). lwd = 0.3 therefore gives a 0.225 pt line, BELOW even
## the 0.3 pt journal minimum -- an easy and silent violation.
##
## The two T&F documents disagree (0.3 pt journals, 0.5 pt CRC books); we design
## to 0.5 pt so both are satisfied, with headroom so the figure survives a
## publisher reduction of ~20%.
PT  <- function(pt) pt / 0.75          # points -> R lwd
LWD_CURVE <- PT(0.7)                   # 0.70 pt: matches measured JASA ROC curves
LWD_AXIS  <- PT(0.6)                   # 0.60 pt
LWD_REF   <- PT(0.5)                   # 0.50 pt: the floor itself

## --- type ------------------------------------------------------------------
## Helvetica: a standard PostScript face T&F name explicitly, and empirically
## what JASA figures use (all 8 sampled papers pair a serif body with
## Helvetica/Arial figures). VERIFIED to resolve on this machine and to embed
## as "Helvetica" -- a font family that does not resolve is silently replaced
## by BitstreamVeraSans with no error and no warning, so this is checked in
## analysis/R/04_check.R rather than assumed.
FONT_FAMILY <- "Helvetica"

## ASA Style Guide: lettering "no smaller than 8 points" at final printed size.
## Every cex multiplier below is >= 1, so 8 pt is the floor for ALL text.
PS_POINTSIZE <- 8
CEX_AXIS   <- 1.00                     # -> 8.0 pt tick labels
CEX_LAB    <- 1.05                     # -> 8.4 pt axis titles
CEX_STRIP  <- 1.00                     # -> 8.0 pt panel identifiers
CEX_LEGEND <- 1.00

## --- canvas ----------------------------------------------------------------
## Measured from typeset JASA articles: full text width 516 pt = 7.17 in.
WIDTH_FULL   <- 7.17
WIDTH_COLUMN <- 3.50

## --- devices ---------------------------------------------------------------

## WHY cairo_ps() AND NOT postscript(), AND WHY NOT embedFonts().
## Measured on this machine, inspecting the EPS natively (not via a ghostscript
## conversion, which fabricates evidence in both directions):
##
##   device                         fonts embedded   gs must load externally
##   postscript(family="Times")     NO               NimbusRoman x4, Symbol
##   postscript() + embedFonts()    partial          NimbusRoman, Symbol
##   cairo_ps()                     YES (Type 42)    nothing
##
## postscript() writes only "%%DocumentNeededResources: font Times-Roman" and
## no font programs -- exactly the case T&F warn about ("any fonts that are not
## embedded will be replaced by Courier"). grDevices::embedFonts() does NOT
## repair it: in four configurations it still failed to embed Symbol (which is
## what plotmath uses for Greek), and its eps2write output additionally
## contains a Type 3 font and a raster image while DELETING the
## %%DocumentNeededResources line -- so a DSC-based check then passes
## vacuously. cairo_ps() embeds a real TrueType subset and needs nothing.
##
## postscript() has a second, worse property: it silently DROPS semi-transparent
## objects (one easily-missed warning, and the file gets smaller). We use no
## alpha anywhere -- T&F forbid transparency and EPS cannot carry it -- but
## cairo_ps rasterises rather than deletes, which is the safer failure mode.
## fallback_resolution = 1200 matches the line-art raster requirement in case
## any transparent object ever slips in.

#' Draw a figure to the submission device and to a screen preview.
#'
#' Produces:
#'   <stem>.eps  -- THE submission file (vector, live text, fonts embedded)
#'   <stem>.png  -- a screen preview for us; NEVER submitted
#'
#' `draw` must be a FUNCTION of no arguments; it is called once per device.
#' (It cannot be an unevaluated expression: an R promise evaluates only once,
#' so the second device would receive an empty page.)
#'
#' `width`/`height` are inches AT FINAL PRINTED SIZE. Do not scale in LaTeX --
#' scaling a figure down shrinks its effective line weights and font sizes
#' proportionally, which is how a compliant figure silently drops below the
#' 0.3 pt line and 8 pt type minima.
save_figure <- function(stem, width, height, draw, dir = FIG_DIR, preview = TRUE,
                        check = TRUE) {
  stopifnot(is.function(draw))
  eps <- file.path(dir, paste0(stem, ".eps"))

  grDevices::cairo_ps(
    filename            = eps,
    width               = width,
    height              = height,
    pointsize           = PS_POINTSIZE,
    family              = FONT_FAMILY,
    onefile             = FALSE,   # EPS must hold exactly one page
    fallback_resolution = 1200,    # = the line-art raster rule, if ever needed
    bg                  = "white"
  )
  tryCatch(draw(), finally = grDevices::dev.off())

  if (preview) {
    png <- file.path(dir, paste0(stem, ".png"))
    grDevices::png(png, width = width, height = height, units = "in",
                   res = 400, pointsize = PS_POINTSIZE, family = FONT_FAMILY,
                   type = "cairo", bg = "white")
    tryCatch(draw(), finally = grDevices::dev.off())
  }

  message("wrote ", eps, if (preview) "  (+ preview .png)" else "")

  ## A figure that fails the artwork spec must not reach the submission
  ## silently. The gate runs on the PRODUCED FILE, because several failure
  ## modes (an unresolved font family, a stray transparent object, a line
  ## thinned below the minimum by a panel layout) produce no R warning at all.
  if (check && exists("check_figure", mode = "function")) {
    rep <- check_figure(eps, width_in = width, height_in = height,
                        min_line_pt = 0.5, scale = 1, min_text_pt = 8)
    if (!isTRUE(rep$ok)) warning("figure ", basename(eps),
                             " does NOT meet the artwork spec (see report above)",
                             call. = FALSE)
  }
  invisible(eps)
}

## --- shared panel furniture ------------------------------------------------

#' An L-shaped axis frame: left and bottom spines only, ticks outward.
#'
#' Both a full box and an L-shape are published JASA practice; the L-shape is
#' what this group's own papers use (Tian et al.'s Comp-gLASSO ROC panels and
#' Tian-Jiang-Jiang's reference-invariance Figure 2), so it is the house style
#' the PI will expect.
panel_axes <- function(xat, yat, xlab_show = TRUE, ylab_show = TRUE,
                       xfmt = NULL, yfmt = NULL) {
  if (is.null(xfmt)) xfmt <- format(xat, trim = TRUE)
  if (is.null(yfmt)) yfmt <- format(yat, trim = TRUE)
  axis(1, at = xat, labels = if (xlab_show) xfmt else FALSE,
       lwd = 0, lwd.ticks = LWD_AXIS, tcl = -0.20,
       cex.axis = CEX_AXIS, mgp = c(3, 0.25, 0))
  axis(2, at = yat, labels = if (ylab_show) yfmt else FALSE,
       lwd = 0, lwd.ticks = LWD_AXIS, tcl = -0.20, las = 1,
       cex.axis = CEX_AXIS, mgp = c(3, 0.45, 0))
  # the L: draw the two spines explicitly so there is no top/right border
  usr <- par("usr")
  segments(usr[1], usr[3], usr[2], usr[3], lwd = LWD_AXIS, xpd = NA)
  segments(usr[1], usr[3], usr[1], usr[4], lwd = LWD_AXIS, xpd = NA)
}

#' A one-line panel identifier drawn in the top margin.
#'
#' This is a panel LABEL (which JASA figures carry), not a figure title (which
#' T&F forbid inside the file).
panel_strip <- function(txt) {
  mtext(txt, side = 3, line = 0.15, cex = par("cex") * CEX_STRIP, adj = 0.5)
}

#' WCAG relative luminance of a colour, for contrast decisions.
relative_luminance <- function(hex) {
  v <- grDevices::col2rgb(hex)[, 1] / 255
  lin <- ifelse(v <= 0.03928, v / 12.92, ((v + 0.055) / 1.055)^2.4)
  sum(c(0.2126, 0.7152, 0.0722) * lin)
}

#' Draw a mean +/- se error bar that stays visible on ANY bar fill.
#'
#' An error bar on a filled bar spans two different backgrounds: below the bar
#' top it lies on the fill, above it lies on white paper. A single colour is
#' therefore guaranteed to disappear on one side -- black error bars vanish
#' inside a dark bar, white ones vanish above it. We split the bar at the mean
#' and draw each half in the colour that contrasts with what is behind it.
#'
#' This keeps the focal series black (the only achromatic, colour-vision-safe,
#' maximum-greyscale-ink choice) without hiding its uncertainty.
error_bar <- function(x, mean, se, fill, half_width) {
  if (!is.finite(se) || se <= 0) return(invisible(NULL))
  lo <- max(0, mean - se); hi <- mean + se
  inside <- if (relative_luminance(fill) < 0.35) "white" else "black"

  ## the segment lying ON the bar
  if (lo < mean) {
    segments(x, lo, x, mean, lwd = LWD_AXIS, col = inside)
    segments(x - half_width, lo, x + half_width, lo, lwd = LWD_AXIS, col = inside)
  }
  ## the segment lying on white paper
  segments(x, mean, x, hi, lwd = LWD_AXIS, col = "black")
  segments(x - half_width, hi, x + half_width, hi, lwd = LWD_AXIS, col = "black")
  invisible(NULL)
}

#' The bottom strip of a multi-panel figure: shared x-axis title + one legend.
#'
#' Both live in a dedicated layout cell rather than in the outer margin, so the
#' axis title sits directly beneath the panels and above the legend. Placing it
#' with mtext() in the outer margin puts it BELOW the legend, and placing it at
#' a small negative line collides with the bottom row's tick labels.
#'
#' One legend for the whole figure, unboxed, following the external-legend
#' precedent in Cao, Lin & Li (JASA) Figure 1 -- rather than repeating it in
#' all sixteen panels.
legend_strip <- function(methods, type = c("line", "fill"), xlab = NULL) {
  type <- match.arg(type)
  labs <- METHOD_LABEL[methods]
  par(mar = c(0, 0, 0, 0))
  plot.new()
  plot.window(xlim = c(0, 1), ylim = c(0, 1))
  if (!is.null(xlab)) {
    text(0.5, 0.88, labels = xlab, adj = c(0.5, 1),
         cex = par("cex") * CEX_LAB, xpd = NA)
  }
  yleg <- if (is.null(xlab)) 0.6 else 0.34
  if (type == "line") {
    legend(0.5, yleg, xjust = 0.5, yjust = 0.5, legend = labs,
           col = METHOD_COL[methods], lty = METHOD_LTY[methods],
           lwd = LWD_CURVE, bty = "n", horiz = TRUE, cex = CEX_LEGEND,
           seg.len = 2.4, xpd = NA)
  } else {
    legend(0.5, yleg, xjust = 0.5, yjust = 0.5, legend = labs,
           fill = METHOD_COL[methods], border = "black", bty = "n",
           horiz = TRUE, cex = CEX_LEGEND, xpd = NA)
  }
}

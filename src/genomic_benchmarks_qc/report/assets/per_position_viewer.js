/* Per-position nucleotide viewer.
 *
 * Draws every nucleotide panel plus the coverage strip into ONE canvas, so a
 * single x-transform keeps the panels aligned and makes zoom linked by
 * construction. Reads its data from a <script type="application/json"> sibling,
 * never from a JS string literal, so a class label containing markup cannot
 * break out of the script.
 *
 * The rendering deliberately reproduces plot_per_base_sequence_comparison so
 * this panel does not look like a foreign object next to the static plots in
 * the same report. Every constant in LAYOUT and STYLE below was measured from
 * that function's output rather than chosen by eye - see the comments for what
 * each one came from. There is one intentional deviation, marked at the point
 * it happens: a flag band never renders thinner than 2px.
 */
(function (global) {
  'use strict';

  /* Geometry in the matplotlib figure's own pixels, measured off its output
   * (3346x2624 for four bases, from a 12x10in figure at dpi 300 after savefig's
   * tight crop). Everything on screen is this, scaled. */
  var LAYOUT = {
    bodyX0: 211,                 // left spine
    bodyX1: 3003,                // right spine
    panelH: 442.5,               // one nucleotide panel
    gapH: 80,                    // between panels (hspace)
    topPad: 29,                  // above the first panel
    bottomPad: 285,              // x labels, axis title and legend
    covRatio: 0.5,               // height_ratios=[1]*n + [0.5]
  };

  /* Where the report's other full-width plots put their axes body, as fractions
   * of the content column: they render at 100% and are 3032px wide with spines
   * at 211 and 3003. The per-position body is sized to land on exactly these,
   * so plot bodies line up down the page. Its wider margins then fall outside
   * that span - the y label hangs to the left, uncounted, which is where the
   * other plots' y labels sit too, because matching the body width also matches
   * the pixel scale (3346 / 1.10356 = 3032). */
  var ALIGN = { bodyX0: 211 / 3032, bodyX1: 3003 / 3032 };

  var DPI_PT = 300 / 72;         // figure pixels per point

  var STYLE = {
    yMax: 1.1,                   // ylim=(0, 1.1) - the shared coverage limit
    xMarginFrac: 0.05,           // matplotlib's default x margin
    font: '"DejaVu Sans", Verdana, Arial, sans-serif',
    labelPt: 14,                 // ylabel / xlabel fontsize=14
    tickPt: 12,                  // tick_params labelsize=12
    legendPt: 10,                // rcParams default; see note in the report
    seriesPt: 1.5,               // rcParams lines.linewidth
    coveragePt: 2,               // plot(..., linewidth=2)
    spinePt: 0.8,                // rcParams axes.linewidth
    tickLenPt: 3.5,              // rcParams xtick.major.size
    tickPadPt: 3.5,              // rcParams xtick.major.pad
    labelPadPt: 4,               // Axis.labelpad - the gap after the tick labels
    seriesAlpha: 0.7,            // plot(..., alpha=0.7)
    flagAlpha: 0.5,              // axvspan(..., alpha=0.5)
    coverageFill: '#d3d3d3',     // 'lightgray'
    coverageAlpha: 0.5,
    /* The static plot wraps this onto three lines. Two keeps the right margin
     * inside the 8% the report's padding allows, which is what lets the body
     * align with the plots above instead of being squeezed narrower. */
    covLabel: ['Proportion of each class', 'reaching position'],
    spine: '#000000',
    text: '#000000',
    legendEdge: '#cccccc',       // rcParams legend.edgecolor '0.8'
    legendFace: '#ffffff',
    paper: '#ffffff',
    unknown: '#9aa4b2',
    minFlagPx: 2,                // the one deviation from the static plot
  };

  var MIN_SPAN = 8;
  var FLAG_WINDOW = 25;

  function niceTicks(x0, x1, target) {
    var span = x1 - x0;
    if (!(span > 0)) return [];
    var rough = span / target;
    var mag = Math.pow(10, Math.floor(Math.log10(rough)));
    var norm = rough / mag;
    var step = (norm <= 1 ? 1 : norm <= 2 ? 2 : norm <= 5 ? 5 : 10) * mag;
    var ticks = [];
    for (var t = Math.ceil(x0 / step) * step; t <= x1; t += step) {
      var v = Math.round(t);
      if (v >= 1 && ticks[ticks.length - 1] !== v) ticks.push(v);
    }
    return ticks;
  }

  /* The tick set plot_per_base_sequence_comparison builds: position 1, then
   * every tenth of the window. Reproduced so the view at rest is identical. */
  function matplotlibTicks(end) {
    var step = Math.max(1, Math.floor(end / 10));
    var ticks = [1];
    for (var t = step; t <= end; t += step) ticks.push(t);
    return ticks;
  }

  function fmt(v) {
    return v == null || isNaN(v) ? '—' : v.toFixed(3);
  }

  function createViewer(root) {
    var dataEl = document.getElementById(root.getAttribute('data-payload'));
    if (!dataEl) return null;
    var data;
    try {
      data = JSON.parse(dataEl.textContent);
    } catch (err) {
      return null; // leave the static fallback in place
    }

    var nts = data.nucleotides || [];
    if (!nts.length || !data.endPosition) return null;

    // Every position in the payload was compared: the figure draws the window
    // the report scored and stops there, so nothing on screen needs explaining
    // away as un-scored.
    var end = data.endPosition;
    var canvas = root.querySelector('.ppv-canvas');
    var tooltip = root.querySelector('.ppv-tooltip');
    var readout = root.querySelector('.ppv-readout');
    var ctx = canvas.getContext('2d');

    // Findings worth navigating to. Unknown should not occur inside the drawn
    // window - it would mean the results table is missing a check - and it is a
    // statement about the cohort rather than a difference between the classes,
    // so where it does turn up it is drawn but left out of the flag count and
    // the jump order.
    var flagList = [];
    Object.keys(data.flags || {}).forEach(function (nt) {
      Object.keys(data.flags[nt]).forEach(function (pos) {
        var flag = data.flags[nt][pos];
        if (flag === 'Unknown') return;
        flagList.push({ pos: +pos, nt: nt, flag: flag });
      });
    });
    flagList.sort(function (a, b) { return a.pos - b.pos; });

    // Merge any Unknown positions into spans once rather than stroking adjacent
    // one-pixel bands on every redraw.
    var unknownRuns = {};
    Object.keys(data.flags || {}).forEach(function (nt) {
      var positions = Object.keys(data.flags[nt])
        .filter(function (p) { return data.flags[nt][p] === 'Unknown'; })
        .map(Number).sort(function (a, b) { return a - b; });
      var runs = [];
      positions.forEach(function (p) {
        var last = runs[runs.length - 1];
        if (last && p === last[1] + 1) last[1] = p;
        else runs.push([p, p]);
      });
      if (runs.length) unknownRuns[nt] = runs;
    });

    /* Coverage arrives per class. The lower of the two curves is what a
     * comparison at a position actually has to work with - a position is scored
     * only where both classes clear the floor - so it is the one the panel
     * fills under. */
    var coverage = data.coverage || [[], []];
    var bindingCoverage = (coverage[0] || []).map(function (value, index) {
      var other = (coverage[1] || [])[index];
      return other == null ? value : Math.min(value, other);
    });

    /* The flag the report gives this base at this position. Only flagged
     * positions travel in the payload, so a position with no entry of its own is
     * one that passed - the hover card names it rather than leaving a blank, so
     * every row there says what the check decided. */
    function flagAt(nt, pos) {
      return (data.flags[nt] || {})[String(pos)] || 'Pass';
    }

    var fullView = {
      x0: 1 - STYLE.xMarginFrac * (end - 1),
      x1: end + STYLE.xMarginFrac * (end - 1),
    };
    var view = { x0: fullView.x0, x1: fullView.x1 };
    var hover = null, hoverPanel = -1, drag = null;
    var colW = 0;      // the content column the plot aligns to
    var G = null;      // layout for that column, rebuilt on resize

    /* ---- geometry ---- */

    /* Build the layout for a content column of `colW`, using `c` only to measure
     * text. The axes body is sized to ALIGN, then the margins are whatever the
     * labels need at that scale, so the canvas ends up wider than the column and
     * is offset to bring the body onto the aligned span. */
    function geom(colW, c) {
      var bodyPx = (ALIGN.bodyX1 - ALIGN.bodyX0) * colW;
      var S = bodyPx / (LAYOUT.bodyX1 - LAYOUT.bodyX0);   // screen px per figure px
      var g = {
        S: S,
        s: function (pt) { return pt * DPI_PT * S; },      // points -> px
        left: LAYOUT.bodyX0 * S,
        plotW: bodyPx,
        panelH: LAYOUT.panelH * S,
        gapH: LAYOUT.gapH * S,
        covH: LAYOUT.panelH * LAYOUT.covRatio * S,
      };
      g.top = LAYOUT.topPad * S;
      g.panelTop = function (i) { return g.top + i * (g.panelH + g.gapH); };
      g.covTop = g.panelTop(nts.length - 1) + g.panelH + g.gapH;
      g.H = g.covTop + g.covH + LAYOUT.bottomPad * S;
      // Right margin: tick marks, tick labels and the rotated coverage label.
      g.rightMargin = g.s(STYLE.tickLenPt) + g.s(STYLE.tickPadPt)
        + tickLabelWidth(c, g) + g.s(STYLE.labelPadPt)
        + STYLE.covLabel.length * g.s(STYLE.labelPt) * 1.2;
      g.W = g.left + g.plotW + g.rightMargin;
      // Shift so the body starts on the aligned span rather than at the canvas
      // edge. Zero whenever the figures share a left margin, which they do.
      g.offsetLeft = ALIGN.bodyX0 * colW - g.left;
      g.xToPx = function (pos) {
        return g.left + (pos - view.x0) / (view.x1 - view.x0) * g.plotW;
      };
      g.pxToX = function (px) {
        return view.x0 + (px - g.left) / g.plotW * (view.x1 - view.x0);
      };
      g.yToPx = function (v, top, h) { return top + h - (v / STYLE.yMax) * h; };
      return g;
    }

    /* Column width that renders the figure at its own pixel scale, for export. */
    function exportColWidth() {
      return (LAYOUT.bodyX1 - LAYOUT.bodyX0) / (ALIGN.bodyX1 - ALIGN.bodyX0);
    }

    function span() { return view.x1 - view.x0; }
    function atView(v) {
      return Math.abs(view.x0 - v.x0) < 1e-6 && Math.abs(view.x1 - v.x1) < 1e-6;
    }
    function atRest() { return atView(fullView); }

    function clampView(x0, x1) {
      if (x1 - x0 < MIN_SPAN) {
        var mid = (x0 + x1) / 2;
        x0 = mid - MIN_SPAN / 2;
        x1 = mid + MIN_SPAN / 2;
      }
      if (x1 - x0 >= fullView.x1 - fullView.x0) {
        x0 = fullView.x0; x1 = fullView.x1;
      }
      if (x0 < fullView.x0) { x1 += fullView.x0 - x0; x0 = fullView.x0; }
      if (x1 > fullView.x1) { x0 -= x1 - fullView.x1; x1 = fullView.x1; }
      view.x0 = Math.max(fullView.x0, x0);
      view.x1 = Math.min(fullView.x1, x1);
      updateReadout();
    }

    function resize() {
      colW = root.querySelector('.ppv-plot').clientWidth;
      if (!colW) return;                     // hidden, e.g. a closed <details>
      G = geom(colW, ctx);
      var dpr = global.devicePixelRatio || 1;
      canvas.width = Math.round(G.W * dpr);
      canvas.height = Math.round(G.H * dpr);
      canvas.style.width = G.W + 'px';
      canvas.style.height = G.H + 'px';
      // the canvas is wider than its column; pull it left so the axes body
      // lands on the aligned span and the y label hangs outside it
      canvas.style.marginLeft = G.offsetLeft + 'px';
      ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
      draw();
    }

    /* ---- drawing primitives that mirror the matplotlib axes ---- */

    function setFont(c, g, pt) {
      c.font = g.s(pt) + 'px ' + STYLE.font;
    }

    function drawSpines(c, g, top, h) {
      c.strokeStyle = STYLE.spine;
      c.lineWidth = g.s(STYLE.spinePt);
      c.strokeRect(
        Math.round(g.left) + 0.5, Math.round(top) + 0.5,
        Math.round(g.plotW), Math.round(h)
      );
    }

    var Y_TICKS = ['0', '0.5', '1'];

    /* Width of the widest y tick label, which is what the axis label has to
     * clear. Measured rather than assumed, because it depends on the font the
     * browser actually resolved. */
    function tickLabelWidth(c, g) {
      setFont(c, g, STYLE.tickPt);
      return Y_TICKS.reduce(function (w, t) {
        return Math.max(w, c.measureText(t).width);
      }, 0);
    }

    /* Centre line for a rotated axis label, laid out the way matplotlib stacks
     * the left margin: tick length, tick pad, the tick labels, labelpad, then
     * the label's own height. Verified against the static PNG, where
     * 'Frequency' sits at x 37..91 against a spine at 212. */
    function axisLabelX(c, g, right, nLines) {
      var block = nLines * g.s(STYLE.labelPt) * 1.2;
      var out = g.s(STYLE.tickLenPt) + g.s(STYLE.tickPadPt) + tickLabelWidth(c, g)
        + g.s(STYLE.labelPadPt) + block / 2;
      // The margins come from a figure measured with DejaVu Sans. A browser that
      // falls back to a wider face makes the tick labels wider, so clamp to keep
      // the label on the canvas even if that eats into the labelpad.
      if (right) return Math.min(g.left + g.plotW + out, g.left + g.plotW + g.rightMargin - block / 2);
      return Math.max(g.left - out, block / 2);
    }

    function drawYAxis(c, g, top, h, right) {
      // yticks [0, 0.5, 1] with outward ticks, as set on the shared axis.
      var len = g.s(STYLE.tickLenPt);
      var pad = g.s(STYLE.tickPadPt);
      c.strokeStyle = STYLE.spine;
      c.lineWidth = g.s(STYLE.spinePt);
      c.fillStyle = STYLE.text;
      setFont(c, g, STYLE.tickPt);
      c.textBaseline = 'middle';
      c.textAlign = right ? 'left' : 'right';
      [0, 0.5, 1].forEach(function (v, i) {
        var y = Math.round(g.yToPx(v, top, h)) + 0.5;
        var x = right ? g.left + g.plotW : g.left;
        c.beginPath();
        c.moveTo(x, y);
        c.lineTo(x + (right ? len : -len), y);
        c.stroke();
        c.fillText(Y_TICKS[i], x + (right ? len + pad : -len - pad), y);
      });
    }

    function drawXAxis(c, g, top, h, ticks) {
      var len = g.s(STYLE.tickLenPt);
      var y = top + h;
      c.strokeStyle = STYLE.spine;
      c.lineWidth = g.s(STYLE.spinePt);
      c.fillStyle = STYLE.text;
      setFont(c, g, STYLE.tickPt);
      c.textAlign = 'center';
      c.textBaseline = 'top';
      ticks.forEach(function (t) {
        var px = g.xToPx(t);
        if (px < g.left - 1 || px > g.left + g.plotW + 1) return;
        px = Math.round(px) + 0.5;
        c.beginPath();
        c.moveTo(px, y);
        c.lineTo(px, y + len);
        c.stroke();
        c.fillText(String(t), px, y + len + g.s(2));
      });
    }

    function drawRotatedLabel(c, g, lines, x, cy, pt) {
      setFont(c, g, pt);
      c.fillStyle = STYLE.text;
      c.save();
      c.translate(x, cy);
      c.rotate(-Math.PI / 2);
      c.textAlign = 'center';
      c.textBaseline = 'middle';
      var lh = g.s(pt) * 1.2;
      var offset = -(lines.length - 1) / 2 * lh;
      lines.forEach(function (line, i) {
        c.fillText(line, 0, offset + i * lh);
      });
      c.restore();
    }

    function drawSeries(c, g, values, color, top, h, pt) {
      var ppp = g.plotW / span();
      var y = function (v) { return g.yToPx(v, top, h); };
      c.save();
      c.beginPath();
      c.rect(g.left, top, g.plotW, h);   // stay inside the axes, as mpl clips
      c.clip();
      c.strokeStyle = color;
      c.lineWidth = g.s(pt == null ? STYLE.seriesPt : pt);
      c.lineJoin = 'round';
      c.globalAlpha = STYLE.seriesAlpha;
      c.beginPath();
      if (ppp >= 1) {
        var from = Math.max(1, Math.floor(view.x0));
        var to = Math.min(end, Math.ceil(view.x1));
        var started = false;
        for (var p = from; p <= to; p++) {
          var v = values[p - 1];
          if (v == null) continue;
          var px = g.xToPx(p);
          if (!started) { c.moveTo(px, y(v)); started = true; }
          else c.lineTo(px, y(v));
        }
      } else {
        // More positions than pixels: a min/max envelope per pixel column,
        // which is both faster and more honest than overplotting.
        for (var col = 0; col <= g.plotW; col++) {
          var a = Math.max(1, Math.ceil(g.pxToX(g.left + col) - 0.5));
          var b = Math.min(end, Math.floor(g.pxToX(g.left + col + 1) + 0.5));
          var lo = Infinity, hi = -Infinity;
          for (var q = a; q <= b; q++) {
            var vv = values[q - 1];
            if (vv == null) continue;
            if (vv < lo) lo = vv;
            if (vv > hi) hi = vv;
          }
          if (lo === Infinity) continue;
          var x = g.left + col + 0.5;
          c.moveTo(x, y(lo));
          c.lineTo(x, y(hi));
        }
      }
      c.stroke();
      c.restore();
    }

    function drawFlags(c, g, nt, top, h) {
      var flags = (data.flags || {})[nt];
      if (!flags) return;
      var ppp = g.plotW / span();
      c.save();
      c.beginPath();
      c.rect(g.left, top, g.plotW, h);
      c.clip();

      // Any Unknown span sits underneath as a flat wash: nothing was concluded
      // there, so it reads as absence of data rather than as a finding.
      (unknownRuns[nt] || []).forEach(function (run) {
        if (run[1] < view.x0 - 1 || run[0] > view.x1 + 1) return;
        var l = g.xToPx(run[0] - 0.5), r = g.xToPx(run[1] + 0.5);
        c.globalAlpha = 0.3;
        c.fillStyle = data.flagColors.Unknown || STYLE.unknown;
        c.fillRect(l, top, Math.max(STYLE.minFlagPx, r - l), h);
      });

      Object.keys(flags).forEach(function (posKey) {
        var flag = flags[posKey];
        if (flag === 'Unknown') return;
        var pos = +posKey;
        if (pos < view.x0 - 1 || pos > view.x1 + 1) return;
        // axvspan(pos - 0.5, pos + 0.5) is one data unit wide, which is
        // sub-pixel on a wide window. This is the one place the drawing departs
        // from the static plot: a band never renders thinner than 2px, because
        // an invisible flag is the problem this view exists to fix.
        var w = Math.max(STYLE.minFlagPx, ppp);
        c.globalAlpha = STYLE.flagAlpha;
        c.fillStyle = data.flagColors[flag] || data.flagColors.Fail;
        c.fillRect(g.xToPx(pos) - w / 2, top, w, h);
      });
      c.restore();
    }

    /* The cohort floor: the sequences that have to reach a position before it
     * may be flagged, as a share of a class. Drawn across the panel because it is
     * the rule that ends the figure: the lower curve reaches it at the right-hand
     * edge, which is what makes that edge read as a limit rather than as the
     * sequences running out.
     *
     * A bare line, with no caption of its own - it used to carry one in italics
     * inside the panel, which put a sentence on top of the curves it explained.
     * The section's ? explanation says what the line is, once. */
    function drawCoverageFloor(c, g, top, h) {
      var floor = data.coverageFloor;
      if (!(floor > 0) || floor > STYLE.yMax) return;
      var y = Math.round(g.yToPx(floor, top, h)) + 0.5;
      c.save();
      c.strokeStyle = 'rgba(0,0,0,0.45)';
      c.lineWidth = 1;
      c.setLineDash([4, 3]);
      c.beginPath();
      c.moveTo(g.left, y);
      c.lineTo(g.left + g.plotW, y);
      c.stroke();
      c.restore();
    }

    /* The bottom panel. One pooled grey curve only said that sequences end
     * somewhere along here; per class, against the floor, it says which class
     * runs out first and why the comparison stops where it does. The grey fill
     * is under the lower of the two curves - the cohort a comparison at that
     * position actually has - and the curves are the two classes in their own
     * colours, so the legend below covers this panel too. */
    function drawCoverage(c, g, ticks) {
      var top = g.covTop, h = g.covH;
      c.save();
      c.beginPath();
      c.rect(g.left, top, g.plotW, h);
      c.clip();
      var from = Math.max(1, Math.floor(view.x0));
      var to = Math.min(end, Math.ceil(view.x1));
      var step = Math.max(1, Math.floor((to - from) / Math.max(1, g.plotW)));
      c.globalAlpha = STYLE.coverageAlpha;
      c.fillStyle = STYLE.coverageFill;
      c.beginPath();
      c.moveTo(g.xToPx(from), g.yToPx(0, top, h));
      for (var p = from; p <= to; p += step) {
        var v = bindingCoverage[p - 1];
        if (v != null) c.lineTo(g.xToPx(p), g.yToPx(v, top, h));
      }
      c.lineTo(g.xToPx(Math.min(to, end)), g.yToPx(0, top, h));
      c.closePath();
      c.fill();
      c.globalAlpha = 1;
      c.restore();

      drawCoverageFloor(c, g, top, h);
      // Second class first, as in the panels above, so the first label's curve
      // is the one on top wherever they cross.
      drawSeries(c, g, coverage[1] || [], data.colors[1], top, h, STYLE.coveragePt);
      drawSeries(c, g, coverage[0] || [], data.colors[0], top, h, STYLE.coveragePt);

      drawSpines(c, g, top, h);
      drawYAxis(c, g, top, h, true);      // ticks and label on the right
      drawXAxis(c, g, top, h, ticks);
      drawRotatedLabel(c, g, STYLE.covLabel,
        axisLabelX(c, g, true, STYLE.covLabel.length), top + h / 2, STYLE.labelPt);
    }

    function drawLegend(c, g) {
      var pt = STYLE.legendPt;
      setFont(c, g, pt);
      var pad = g.s(4), handle = g.s(pt * 2), gapHL = g.s(4), between = g.s(12);
      var entries = [
        { color: data.colors[0], label: data.labels[0] },
        { color: data.colors[1], label: data.labels[1] },
      ];
      var w = pad * 2;
      entries.forEach(function (e, i) {
        w += handle + gapHL + c.measureText(e.label).width + (i ? between : 0);
      });
      var h = g.s(pt) * 1.6 + pad * 2;
      var x = g.left + g.plotW / 2 - w / 2;
      /* Under the x axis title, in the space the figure's bottom margin leaves
       * below it, and clamped so it cannot run off the bottom of the canvas -
       * which is where it used to end up, laid out against the canvas width
       * instead of the figure scale. */
      var below = g.covTop + g.covH + g.s(STYLE.tickLenPt) + g.s(STYLE.tickPt) * 1.6
        + g.s(STYLE.labelPt) * 1.2 + g.s(4);
      var y = Math.min(below, g.H - h - g.s(2));
      var r = g.s(3);

      c.globalAlpha = 0.8;                 // rcParams legend.framealpha
      c.fillStyle = STYLE.legendFace;
      c.strokeStyle = STYLE.legendEdge;
      c.lineWidth = g.s(STYLE.spinePt);
      c.beginPath();
      c.moveTo(x + r, y);
      c.arcTo(x + w, y, x + w, y + h, r);
      c.arcTo(x + w, y + h, x, y + h, r);
      c.arcTo(x, y + h, x, y, r);
      c.arcTo(x, y, x + w, y, r);
      c.closePath();
      c.fill();
      c.stroke();
      c.globalAlpha = 1;

      var cx = x + pad, cy = y + h / 2;
      c.textBaseline = 'middle';
      c.textAlign = 'left';
      entries.forEach(function (e, i) {
        if (i) cx += between;
        c.strokeStyle = e.color;
        c.lineWidth = g.s(STYLE.seriesPt);
        c.globalAlpha = STYLE.seriesAlpha;
        c.beginPath();
        c.moveTo(cx, cy);
        c.lineTo(cx + handle, cy);
        c.stroke();
        c.globalAlpha = 1;
        cx += handle + gapHL;
        c.fillStyle = STYLE.text;
        c.fillText(e.label, cx, cy);
        cx += c.measureText(e.label).width;
      });
    }

    function drawHover(c, g) {
      if (hover == null) return;
      var px = Math.round(g.xToPx(hover)) + 0.5;
      c.strokeStyle = 'rgba(0,0,0,0.35)';
      c.lineWidth = 1;
      c.setLineDash([3, 3]);
      c.beginPath();
      c.moveTo(px, g.top);
      c.lineTo(px, g.covTop + g.covH);
      c.stroke();
      c.setLineDash([]);
      function dot(value, top, h, series) {
        if (value == null) return;
        c.beginPath();
        c.arc(g.xToPx(hover), g.yToPx(value, top, h), g.s(2.4), 0, Math.PI * 2);
        c.fillStyle = data.colors[series];
        c.fill();
        c.lineWidth = g.s(1.2);          // surface ring keeps the two dots
        c.strokeStyle = STYLE.paper;     // readable where they overlap
        c.stroke();
      }
      nts.forEach(function (nt, i) {
        var top = g.panelTop(i);
        [0, 1].forEach(function (s) { dot(data.freq[nt][s][hover - 1], top, g.panelH, s); });
      });
      [0, 1].forEach(function (s) {
        dot((coverage[s] || [])[hover - 1], g.covTop, g.covH, s);
      });
    }

    function drawSelection(c, g) {
      if (!drag || drag.panning || drag.currentPx == null) return;
      var a = Math.min(drag.startPx, drag.currentPx);
      var b = Math.max(drag.startPx, drag.currentPx);
      if (b - a < 2) return;
      c.fillStyle = 'rgba(0,61,153,0.10)';
      c.fillRect(a, g.top, b - a, g.covTop + g.covH - g.top);
      c.strokeStyle = 'rgba(0,61,153,0.55)';
      c.lineWidth = 1;
      c.beginPath();
      c.moveTo(a + 0.5, g.top); c.lineTo(a + 0.5, g.covTop + g.covH);
      c.moveTo(b - 0.5, g.top); c.lineTo(b - 0.5, g.covTop + g.covH);
      c.stroke();
    }

    /* One render pass. Used for the live canvas and, at the static plot's own
     * scale and without the interaction overlays, for the PNG export. */
    function render(c, g, opts) {
      var W = g.W;
      var ticks = atRest() ? matplotlibTicks(end) : niceTicks(view.x0, view.x1, 10);
      c.fillStyle = STYLE.paper;
      c.fillRect(0, 0, W, g.H);
      nts.forEach(function (nt, i) {
        var top = g.panelTop(i);
        drawFlags(c, g, nt, top, g.panelH);
        drawSeries(c, g, data.freq[nt][1], data.colors[1], top, g.panelH);
        drawSeries(c, g, data.freq[nt][0], data.colors[0], top, g.panelH);
        drawSpines(c, g, top, g.panelH);
        drawYAxis(c, g, top, g.panelH, false);
        drawXAxis(c, g, top, g.panelH, ticks);
        drawRotatedLabel(c, g, ['Frequency'],
          axisLabelX(c, g, false, 1), top + g.panelH / 2, STYLE.labelPt);
        setFont(c, g, STYLE.labelPt);
        c.fillStyle = STYLE.text;
        c.textAlign = 'center';
        c.textBaseline = 'bottom';
        // text(0.9, 0.8, f'Nucleotide: {nt}', transform=axes)
        c.fillText('Nucleotide: ' + nt,
          g.left + 0.9 * g.plotW, top + g.panelH - 0.8 * g.panelH);
      });
      drawCoverage(c, g, ticks);
      setFont(c, g, STYLE.labelPt);
      c.fillStyle = STYLE.text;
      c.textAlign = 'center';
      c.textBaseline = 'top';
      c.fillText(data.xLabel || 'Position in sequence',
        g.left + g.plotW / 2,
        g.covTop + g.covH + g.s(STYLE.tickLenPt) + g.s(STYLE.tickPt) * 1.6);
      drawLegend(c, g);
      if (opts && opts.overlay) {
        drawHover(c, g);
        drawSelection(c, g);
      }
      return g;
    }

    function draw() {
      if (!G) return;
      render(ctx, G, { overlay: true });
    }

    /* ---- export ---- */

    function visibleRange() {
      return [Math.max(1, Math.round(view.x0 + 0.5)), Math.min(end, Math.round(view.x1 - 0.5))];
    }

    function savePng() {
      // Export at the static plot's own pixel scale so a saved crop drops into
      // a figure beside the report's PNGs without resampling.
      var off = document.createElement('canvas');
      var octx = off.getContext('2d');
      var g = geom(exportColWidth(), octx);
      off.width = Math.round(g.W);
      off.height = Math.round(g.H);
      render(octx, g, { overlay: false });
      var range = visibleRange();
      var name = (data.direction === 'reversed'
        ? 'per_position_reversed_nucleotide_content'
        : 'per_position_nucleotide_content')
        + '_' + range[0] + '-' + range[1] + '.png';
      var link = document.createElement('a');
      link.download = name;
      link.href = off.toDataURL('image/png');
      document.body.appendChild(link);
      link.click();
      document.body.removeChild(link);
      return name;
    }

    /* ---- readout & tooltip ---- */

    function updateReadout() {
      if (!readout) return;
      var r = visibleRange();
      if (atRest()) {
        readout.textContent = 'Showing all ' + end + ' positions';
      } else {
        readout.textContent = 'Positions ' + r[0] + '–' + r[1] + ' of ' + end;
      }
    }

    /* The cohort behind the numbers above it: the bottom panel read at the
     * hovered position, in sequences as well as per cent, because 'a fifth of
     * the class' means something different in 200 sequences and in 2 million.
     *
     * One class per row rather than one run-on line: class labels are often
     * bare numbers ('0' and '1'), and beside a percentage on the same line
     * there is nothing to tell the label from the measurement. */
    function reachBlock() {
      var rows = [0, 1].map(function (i) {
        var frac = (coverage[i] || [])[hover - 1];
        if (frac == null) return null;
        var count = (data.counts || [])[i];
        // The per cent is bracketed after the count it restates - '450 of 600
        // (75%)' - so the row reads as one measurement rather than as two
        // numbers side by side. Unbracketed when there is no count to restate.
        var reach = count == null
          ? '<td>' + Math.round(frac * 100) + '%</td><td></td>'
          : '<td>' + Math.round(frac * count).toLocaleString()
            + ' of ' + count.toLocaleString() + '</td>'
            + '<td>(' + Math.round(frac * 100) + '%)</td>';
        return '<tr><th><span class="ppv-swatch" style="background:' + data.colors[i]
          + '"></span>' + data.labels[i] + '</th>' + reach + '</tr>';
      }).filter(Boolean);
      return rows.length
        ? '<div class="ppv-tt-reach"><div class="ppv-tt-reach-head">'
          + 'Sequences reaching this position</div><table><tbody>'
          + rows.join('') + '</tbody></table></div>'
        : '';
    }

    function showTooltip(clientX, clientY) {
      if (hover == null) { tooltip.hidden = true; return; }
      /* Every base gets a chip, Pass included. The panels only shade what is
       * wrong, so without it the card would name a flag on some rows and say
       * nothing on the others, and silence there reads as missing data rather
       * than as a position that was checked and passed. */
      var rows = nts.map(function (nt) {
        var v0 = data.freq[nt][0][hover - 1];
        var v1 = data.freq[nt][1][hover - 1];
        var fl = flagAt(nt, hover);
        var chip = '<span class="ppv-chip" style="background:'
          + data.flagColors[fl] + '">' + fl + '</span>';
        return '<tr' + (nts[hoverPanel] === nt ? ' class="is-focus"' : '') + '>'
          + '<th>' + nt + '</th>'
          + '<td>' + fmt(v0) + '</td><td>' + fmt(v1) + '</td>'
          + '<td class="ppv-tt-flag">' + chip + '</td></tr>';
      }).join('');
      // The two numbers are the share of the sequences reaching this position
      // that carry this base here. Saying so in the card is the only place a
      // reader is told what the columns hold.
      tooltip.innerHTML =
        '<div class="ppv-tt-head">Position ' + hover + '</div>'
        + '<div class="ppv-tt-sub">Frequency of each base, per class</div>'
        + '<table><thead><tr><th>Base</th>'
        + '<td><span class="ppv-swatch" style="background:' + data.colors[0] + '"></span>' + data.labels[0] + '</td>'
        + '<td><span class="ppv-swatch" style="background:' + data.colors[1] + '"></span>' + data.labels[1] + '</td>'
        + '<td>Flag</td></tr></thead><tbody>' + rows + '</tbody></table>'
        + reachBlock();
      tooltip.hidden = false;

      var box = root.querySelector('.ppv-plot').getBoundingClientRect();
      var tw = tooltip.offsetWidth, th = tooltip.offsetHeight;
      var x = clientX - box.left + 14, y = clientY - box.top + 14;
      if (x + tw > box.width) x = clientX - box.left - tw - 14;
      if (y + th > box.height) y = Math.max(4, box.height - th - 4);
      tooltip.style.left = Math.max(4, x) + 'px';
      tooltip.style.top = y + 'px';
    }

    /* ---- interaction ---- */

    function localPx(evt) { return evt.clientX - canvas.getBoundingClientRect().left; }

    function panelAt(evt) {
      var g = G;
      var y = evt.clientY - canvas.getBoundingClientRect().top;
      for (var i = 0; i < nts.length; i++) {
        var t = g.panelTop(i);
        if (y >= t && y <= t + g.panelH) return i;
      }
      return -1;
    }

    function inPlot(g, px) { return px >= g.left && px <= g.left + g.plotW; }

    canvas.addEventListener('mousemove', function (evt) {
      if (!G) return;
      var g = G, px = localPx(evt);
      if (drag) {
        if (drag.panning) {
          var dx = (px - drag.startPx) / g.plotW * drag.originSpan;
          clampView(drag.originX0 - dx, drag.originX0 - dx + drag.originSpan);
        } else {
          drag.currentPx = Math.min(g.left + g.plotW, Math.max(g.left, px));
        }
      }
      if (inPlot(g, px)) {
        hover = Math.min(end, Math.max(1, Math.round(g.pxToX(px))));
        hoverPanel = panelAt(evt);
        showTooltip(evt.clientX, evt.clientY);
      } else {
        hover = null;
        tooltip.hidden = true;
      }
      draw();
    });

    canvas.addEventListener('mouseleave', function () {
      hover = null;
      tooltip.hidden = true;
      draw();
    });

    canvas.addEventListener('mousedown', function (evt) {
      if (!G) return;
      var g = G, px = localPx(evt);
      if (!inPlot(g, px)) return;
      drag = evt.shiftKey
        ? { startPx: px, panning: true, originX0: view.x0, originSpan: span() }
        : { startPx: px, currentPx: px, panning: false };
      evt.preventDefault();
    });

    global.addEventListener('mouseup', function () {
      if (!drag || !G) return;
      var g = G;
      if (!drag.panning && drag.currentPx != null) {
        var a = Math.min(drag.startPx, drag.currentPx);
        var b = Math.max(drag.startPx, drag.currentPx);
        if (b - a >= 4) clampView(g.pxToX(a), g.pxToX(b));
      }
      drag = null;
      draw();
    });

    canvas.addEventListener('wheel', function (evt) {
      if (!G) return;
      evt.preventDefault();
      var g = G, px = localPx(evt);
      var anchor = g.pxToX(Math.min(g.left + g.plotW, Math.max(g.left, px)));
      var newSpan = span() * (evt.deltaY > 0 ? 1.18 : 1 / 1.18);
      var frac = (anchor - view.x0) / span();
      clampView(anchor - frac * newSpan, anchor + (1 - frac) * newSpan);
      draw();
    }, { passive: false });

    canvas.addEventListener('dblclick', function () {
      clampView(fullView.x0, fullView.x1);
      draw();
    });

    function gotoFlag(dir) {
      if (!flagList.length) return;
      var centre = (view.x0 + view.x1) / 2, target = null, i;
      if (dir > 0) {
        for (i = 0; i < flagList.length; i++) {
          if (flagList[i].pos > centre + 0.5) { target = flagList[i]; break; }
        }
        if (!target) target = flagList[0];
      } else {
        for (i = flagList.length - 1; i >= 0; i--) {
          if (flagList[i].pos < centre - 0.5) { target = flagList[i]; break; }
        }
        if (!target) target = flagList[flagList.length - 1];
      }
      clampView(target.pos - FLAG_WINDOW - 0.5, target.pos + FLAG_WINDOW + 0.5);
      hover = target.pos;
      draw();
    }

    root.querySelectorAll('[data-action]').forEach(function (btn) {
      btn.addEventListener('click', function () {
        var a = btn.getAttribute('data-action');
        if (a === 'reset') clampView(fullView.x0, fullView.x1);
        else if (a === 'next') gotoFlag(1);
        else if (a === 'prev') gotoFlag(-1);
        else if (a === 'save') savePng();
        draw();
      });
    });

    function zoomTo(pos) {
      clampView(pos - FLAG_WINDOW - 0.5, pos + FLAG_WINDOW + 0.5);
      hover = pos;
      draw();
    }

    root.querySelectorAll('[data-goto]').forEach(function (link) {
      link.addEventListener('click', function (evt) {
        evt.preventDefault();
        zoomTo(+link.getAttribute('data-goto'));
      });
    });

    /* The flag list is built here rather than emitted as HTML: every flag is
     * already in the payload because the plot needs it, so rows in the markup
     * would be the same data twice. That also makes the row cap a display
     * choice - nothing is withheld from the page. */
    var SEVERITY = { Fail: 0, Warning: 1 };
    var ROW_CAP = 12;

    function buildFlagList() {
      var host = root.querySelector('.ppv-flags-body');
      var summary = root.querySelector('.ppv-flags-count');
      if (!host) return;

      // `|| 9` would be wrong here: Fail's rank is 0, which is falsy.
      var rank = function (f) {
        return Object.prototype.hasOwnProperty.call(SEVERITY, f.flag) ? SEVERITY[f.flag] : 9;
      };
      var gap = function (f) {
        return Math.abs(data.freq[f.nt][0][f.pos - 1] - data.freq[f.nt][1][f.pos - 1]);
      };
      var rows = flagList.slice().sort(function (a, b) {
        var d = rank(a) - rank(b);
        if (d) return d;
        // Within a severity, the widest gap between the classes first. This is
        // the order the AU-ROC column used to give: for a 0/1 indicator like
        // 'is this base here', AU-ROC is 0.5 + gap / 2, so ranking on the gap
        // and ranking on the score are the same ranking. It orders the rows
        // without being shown - it is the two frequencies, subtracted.
        return gap(b) - gap(a) || a.pos - b.pos;
      });

      if (summary) {
        var counts = {};
        rows.forEach(function (r) { counts[r.flag] = (counts[r.flag] || 0) + 1; });
        var parts = ['Fail', 'Warning'].filter(function (k) { return counts[k]; })
          .map(function (k) { return counts[k] + ' ' + k.toLowerCase(); });
        summary.textContent = rows.length
          ? rows.length + ' flagged position' + (rows.length === 1 ? '' : 's')
            + (parts.length ? ' (' + parts.join(', ') + ')' : '')
          : (data.compared === false ? 'Nothing was compared' : 'No flagged positions');
      }
      if (!rows.length) {
        // An empty list means two opposite things, and reading the wrong one is
        // reading an unscored comparison as a clean one.
        host.innerHTML = data.compared === false
          ? '<p class="ppv-empty">No position here was compared: too few sequences reach any of '
            + 'them. The frequencies are drawn so they can be read by eye, and nothing in the '
            + 'figure is a finding.</p>'
          : '<p class="ppv-empty">Every position passed in this direction.</p>';
        return;
      }

      var shown = Math.min(ROW_CAP, rows.length);
      function cell(r) {
        return '<tr tabindex="0" data-pos="' + r.pos + '">'
          + '<td class="ppv-c-flag"><span class="ppv-sev" style="background:'
          + data.flagColors[r.flag] + '"></span>' + r.flag + '</td>'
          + '<td class="ppv-c-pos">' + r.pos + '</td>'
          + '<td class="ppv-c-base">' + r.nt + '</td>'
          + '<td class="ppv-c-num">' + data.freq[r.nt][0][r.pos - 1].toFixed(3) + '</td>'
          + '<td class="ppv-c-num">' + data.freq[r.nt][1][r.pos - 1].toFixed(3) + '</td>'
          + '<td class="ppv-c-go"><span class="ppv-go">Zoom</span></td></tr>';
      }
      /* One header row, and what the numeric columns hold said above the table
       * rather than in a cell spanning them. The spanning cell was centred over
       * two columns and sat a row above the class names, so it read as a
       * heading for the table and not as a description of those two numbers.
       * The caption is the hover card's subtitle, in the same words: the card
       * and the table describe the same measurement.
       *
       * The columns themselves then only need the class names, because the row
       * already says which base and which position. */
      function swatch(i) {
        return '<span class="ppv-swatch" style="background:' + data.colors[i] + '"></span>';
      }
      var posLabel = data.xLabel || 'Position in sequence';
      function paint(limit) {
        host.innerHTML =
          '<p class="ppv-flags-caption">Frequency of this base at this position, '
          + 'per class</p>'
          + '<table class="ppv-flagtable"><thead>'
          + '<tr><th class="ppv-h-left">Flag</th>'
          + '<th class="ppv-h-left">' + posLabel + '</th>'
          + '<th class="ppv-h-left">Base</th>'
          + '<th>' + swatch(0) + data.labels[0] + '</th>'
          + '<th>' + swatch(1) + data.labels[1] + '</th>'
          + '<th></th></tr></thead><tbody>'
          + rows.slice(0, limit).map(cell).join('')
          + '</tbody></table>'
          + (limit < rows.length
            ? '<button type="button" class="ppv-more">Show all ' + rows.length + ' positions</button>'
            : '');
        host.querySelectorAll('tbody tr').forEach(function (tr) {
          var pos = +tr.getAttribute('data-pos');
          tr.addEventListener('click', function () { zoomTo(pos); });
          tr.addEventListener('keydown', function (evt) {
            if (evt.key === 'Enter' || evt.key === ' ') { evt.preventDefault(); zoomTo(pos); }
          });
        });
        var more = host.querySelector('.ppv-more');
        if (more) more.addEventListener('click', function () { paint(rows.length); });
      }
      paint(shown);
    }
    buildFlagList();

    if (global.ResizeObserver) {
      new global.ResizeObserver(resize).observe(root.querySelector('.ppv-plot'));
    } else {
      global.addEventListener('resize', resize);
    }

    // The markup carries a line of text for readers without JavaScript; drop it
    // now that the figure is drawn.
    var fallback = root.querySelector('.ppv-fallback');
    if (fallback) fallback.remove();
    root.classList.add('is-live');

    resize();
    updateReadout();
    return { draw: draw, resize: resize, savePng: savePng };
  }

  global.initPerPositionViewers = function () {
    document.querySelectorAll('.ppv').forEach(createViewer);
  };
}(window));

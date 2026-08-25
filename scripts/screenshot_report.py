#!/usr/bin/env python3
"""Capture the report screenshot that the documentation landing page uses.

The image was taken by hand until now, which meant it aged without anyone
noticing: the one this replaces still carried a tool description three
revisions old and an input path that no longer existed, while sitting on the
front page as evidence of what the current code produces. Generating it from a
report the build has just built puts it under the same rule as everything else
under docs/reports/ - it cannot show anything the code does not currently do.

The frame is deliberate rather than the top of the page. It scrolls to the
per-position viewer of `hidden-motif`, because that panel is the one thing in a
report that prose alongside it cannot stand in for: 398 positions of flat noise
with a six-position spike in the middle of it. The report's sidebar is
`position: fixed`, so it stays in shot and brings the flag summary with it.

It also records the README's two animations. `--mode scroll` scrolls a report
from the top to the bottom, pausing on each section, which is the one thing
neither a still nor prose gets across: that a report is a long document with a
verdict and a figure on every check in it, not a single number. `--mode panel`
drives the per-position panel instead - drag, tooltip, `Next flag` - because
that panel's whole value is in motion.

None of the three is committed. `build_docs.py` generates all of them into
docs/assets/ from the reports it has just built, and the README points at the
copies the site publishes rather than at paths in the repository, which neither
PyPI nor a branch view resolves. That also keeps them under the same rule as
the reports themselves: an image of the report can never be older than the code
that produced it.

Usage:

    python scripts/screenshot_report.py                    # the landing-page shot
    python scripts/screenshot_report.py --mode scroll      # the README page tour
    python scripts/screenshot_report.py --mode panel       # the README panel demo
    python scripts/screenshot_report.py --report <html> --anchor <element-id>

This needs Playwright and its Chromium build:

    pip install playwright && python -m playwright install chromium

Without them it writes a placeholder and succeeds, the same way
`build_docs.py --skip-reports` writes placeholder reports: working on prose
should not require a 150 MB browser download, and `mkdocs --strict` fails on a
missing image, so something has to be at that path. Pass `--require` to fail
instead - `build_docs.py` does for any build that is not a prose build, so a
placeholder cannot reach the published site by way of a broken CI install.
Both animations require it unconditionally: nothing in the site links them, so
`mkdocs --strict` would not notice a placeholder standing in for one.
"""

import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent

# hidden-motif, and the forward per-position panel within it. See the module
# docstring for why this report and not one with more failures in it.
DEFAULT_REPORT = (ROOT / 'docs' / 'reports' / 'hidden-motif' / 'class'
                  / 'sequence' / '0_vs_1' / 'gb-qc-report.html')
DEFAULT_OUT = ROOT / 'docs' / 'assets' / 'report-screenshot.png'
DEFAULT_ANCHOR = 'per-position-nucleotide-content'

# Beside the landing-page shot, and generated like it: see the module docstring.
DEFAULT_PANEL_OUT = ROOT / 'docs' / 'assets' / 'per-position-demo.webp'
DEFAULT_SCROLL_OUT = ROOT / 'docs' / 'assets' / 'report-scroll.webp'

# The tour scrolls composition-bias, the worst case in the gallery, because it
# is the only report that carries all three flags at once: six checks Fail, one
# Warns, two Pass. Scrolling a cleaner one demonstrates the layout; scrolling
# this one demonstrates what the layout is for. Note the pairing that makes it
# the right choice rather than just the loudest - the two flags a reader most
# needs to tell apart, Warning and Fail, appear on adjacent sections here, on
# duplicates within a class against duplicates between them.
DEFAULT_SCROLL_REPORT = (ROOT / 'docs' / 'reports' / 'composition-bias'
                         / 'class' / 'sequence' / '0_vs_1'
                         / 'gb-qc-report.html')

# Wide enough to clear the report's own 900px breakpoint, below which the
# sidebar stops being a left rail and folds into a row of chips - which would
# lose the flag summary the shot is partly there to show. The scale factor is
# what keeps the plot's hairlines from disappearing when the page renders the
# image at half this width.
VIEWPORT = {'width': 1440, 'height': 940}
SCALE = 2

# Recorded at 1x, and 1280 wide so the panel is laid out the way a desktop
# reader sees it rather than in the report's under-900px layout. No 2x here: an
# animation is watched rather than studied, and the scale factor would treble
# the file for detail nobody pauses on. The height only has to exceed the
# panel's, so that Playwright frames it without scrolling between shots.
PANEL_VIEWPORT = {'width': 1280, 'height': 1000}

# WebP stores a duration per frame, so a pause costs one frame held for most of
# a second instead of a dozen identical ones - which is what keeps a demo this
# long down to about two dozen frames.
PANEL_MOVE_MS = 70          # while something is moving
PANEL_HOLD_MS = 800         # after a view changes, long enough to register
PANEL_READ_MS = 1600        # on the tooltip, which has numbers in it to read
PANEL_END_MS = 1500         # before the loop restarts, so it does not jump-cut

PANEL_DRAG_STEPS = 10

# The zoom the drag lands on, as a half-width in positions either side of the
# finding's centre. Four times the finding's own width makes the zoom a visible
# change from the whole sequence while still separating the flagged columns;
# the floor is what gives a one-position finding a window at all.
PANEL_WINDOW_SPREAD = 4
PANEL_MIN_HALF_WINDOW = 24

# Enough air above the section's card that it does not sit flush against the top
# edge, but less than the gap between cards - any more and the card above starts
# showing at the top of the frame as a stray sliver.
ANCHOR_MARGIN = 14

# The tour is captured at the width the README renders it at, so no resampling
# stands between the report's text and the reader. 1x for the same reason
# --mode panel is.
SCROLL_VIEWPORT = {'width': 1200, 'height': 780}

# The report's own table of contents, which is also the order a reader would go
# through it in. Each entry is one stop.
SCROLL_NAV = '.sidebar-item a[href^="#"]'

# How far the page moves between two frames of a glide, and the cap on how many
# frames one glide can take. Together they mean every jump travels at the same
# speed, except that a very long one speeds up rather than stretching the
# animation - the gap between two sections is not itself interesting.
SCROLL_STEP_PX = 210
SCROLL_MAX_STEPS = 11
# The run back to the top at the end. Fast enough to read as a rewind rather
# than as a second tour, and it is what makes the loop seamless: the last frame
# is one step short of the top, so wrapping round to the first is one more step
# of the same length.
SCROLL_RETURN_PX = 900

SCROLL_GLIDE_MS = 55       # while the page is moving
SCROLL_HOLD_MS = 850       # on a section, long enough to see what it says
SCROLL_TOP_MS = 1400       # at the top, where the loop begins and ends
SCROLL_END_MS = 1700       # on the last section, before the rewind

# The report's own scroll position is not in the frame - headless Chromium draws
# no scrollbar, and the page is 1200px wide because that is how the README shows
# it, not because a scrollbar was wanted there. Without some indicator, a glide
# between two sections that look alike reads as a cut. So one is drawn on.
SCROLL_SETTLE_MS = 40
SCROLL_BAR_WIDTH = 6
SCROLL_BAR_INSET = 9
SCROLL_BAR_MIN = 34
SCROLL_BAR_FILL = (16, 21, 28, 105)

# Chromium's own smooth scrolling, and the report's scroll-behavior, would both
# still be moving when the frame is taken.
PREP_CSS = """
*, *::before, *::after {
    animation: none !important;
    transition: none !important;
}
html, body { scroll-behavior: auto !important; }
"""

# The stops, read off the report's sidebar in document order. Returning the
# document offset rather than the link means the scroll is driven by where each
# section actually is, so a report with a section more or fewer just gets a
# longer or shorter tour.
SCROLL_STOPS = """
selector => {
  const seen = new Set(), out = [];
  for (const link of document.querySelectorAll(selector)) {
    const href = link.getAttribute('href');
    if (!href || href.length < 2 || seen.has(href)) continue;
    const target = document.getElementById(decodeURIComponent(href.slice(1)));
    if (!target) continue;
    seen.add(href);
    out.push({label: (link.textContent || href).trim(),
              y: target.getBoundingClientRect().top + window.scrollY});
  }
  return out;
}
"""

# The per-position panel paints into a <canvas> from a JSON payload after load,
# and an empty canvas screenshots as a white rectangle without erroring - the
# exact failure that would go unnoticed in CI. So wait for the canvas to hold
# more than one colour, which is only true once something has been drawn into
# it. Reports with no such panel have nothing to wait for.
CANVAS_PAINTED = """
() => {
  const canvas = document.querySelector('.ppv-canvas');
  if (!canvas) return true;
  if (!canvas.width || !canvas.height) return false;
  const data = canvas.getContext('2d')
    .getImageData(0, 0, canvas.width, canvas.height).data;
  for (let i = 4; i < data.length; i += 4) {
    if (data[i] !== data[0] || data[i + 1] !== data[1]
        || data[i + 2] !== data[2]) return true;
  }
  return false;
}
"""

# Everything flagged in the panel, straight out of the payload the plot itself
# draws from, so --animate frames whatever the report it is given actually found
# instead of positions written down here.
FLAGGED_POSITIONS = """
root => {
  const el = document.getElementById(root.getAttribute('data-payload'));
  const data = JSON.parse(el.textContent);
  const all = new Set(), fails = new Set();
  Object.keys(data.flags || {}).forEach(function (nt) {
    Object.keys(data.flags[nt]).forEach(function (pos) {
      all.add(+pos);
      if (data.flags[nt][pos] === 'Fail') fails.add(+pos);
    });
  });
  const sorted = s => Array.from(s).sort((a, b) => a - b);
  return {end: data.endPosition, all: sorted(all), fails: sorted(fails)};
}
"""

# Which position the plot thinks the cursor is over. Two of these either side of
# the plot solve its data-to-pixel mapping, which is how the animation aims a
# drag at a position without reproducing the viewer's axis arithmetic here -
# where it could quietly fall out of step with it.
HOVER_POSITION = """
root => {
  const head = root.querySelector('.ppv-tt-head');
  if (!head || root.querySelector('.ppv-tooltip').hidden) return null;
  const n = parseInt(head.textContent.replace(/[^0-9]/g, ''), 10);
  return Number.isFinite(n) ? n : null;
}
"""

# A screenshot does not include the pointer, so an animation driven by one shows
# things happening for no visible reason. This draws a stand-in at the same
# coordinates the real pointer is moved to - the interaction underneath is the
# real thing, only the arrow is painted. Positioned in document coordinates
# rather than fixed, so it stays put against the panel if the page scrolls.
CURSOR_SETUP = """
() => {
  const el = document.createElement('div');
  el.style.cssText = 'position:absolute;z-index:99999;pointer-events:none;'
    + 'width:20px;height:20px;margin:-2px 0 0 -2px;opacity:0';
  el.innerHTML =
    '<svg width="20" height="20" viewBox="0 0 20 20" style="overflow:visible">'
    + '<circle class="ring" cx="6" cy="7" r="11" fill="none" stroke="#1f6fd0"'
    + ' stroke-width="2" opacity="0"/>'
    + '<path d="M1 1 L1 15 L5 11 L7.6 16.8 L10.4 15.5 L7.8 9.9 L13.4 9.6 Z"'
    + ' fill="#10151c" stroke="#ffffff" stroke-width="1.3"'
    + ' stroke-linejoin="round"/></svg>';
  document.body.appendChild(el);
  window.__gbqcCursor = function (x, y, down) {
    el.style.left = (x + window.scrollX) + 'px';
    el.style.top = (y + window.scrollY) + 'px';
    el.style.opacity = '1';
    el.querySelector('.ring').setAttribute('opacity', down ? '0.9' : '0');
  };
}
"""


def write_placeholder(out: Path, reason: str):
    """Write an obviously-fake image, so a real one is never faked silently.

    It has to be a valid PNG at the right path or the strict build fails on the
    missing link, but it must also be unmistakable on sight: if one of these
    ever reaches the site, the mistake should be visible in the first second of
    looking at the page rather than found later.
    """
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    fig = plt.figure(figsize=(VIEWPORT['width'] / 100,
                              VIEWPORT['height'] / 100), dpi=100)
    fig.patch.set_facecolor('#f3b41a')
    fig.text(0.5, 0.56, 'PLACEHOLDER', ha='center', va='center',
             fontsize=52, fontweight='bold', color='#10151c')
    fig.text(0.5, 0.44, f'report screenshot not generated: {reason}',
             ha='center', va='center', fontsize=17, color='#10151c')
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out)
    plt.close(fig)


def capture(report: Path, out: Path, anchor: str):
    """Screenshot one report, framed on `anchor`, with the sidebar in shot."""
    from playwright.sync_api import sync_playwright

    out.parent.mkdir(parents=True, exist_ok=True)
    with sync_playwright() as playwright:
        browser = playwright.chromium.launch()
        page = browser.new_page(viewport=VIEWPORT, device_scale_factor=SCALE)
        page.goto(report.resolve().as_uri())
        page.wait_for_load_state('load')
        page.wait_for_function(CANVAS_PAINTED, timeout=30_000)

        # `instant` matters: the report sets scroll-behavior on some elements,
        # and a smooth scroll would still be moving when the shot is taken.
        found = page.evaluate(
            """([id, margin]) => {
                const target = document.getElementById(id);
                if (!target) return false;
                target.scrollIntoView({block: 'start', behavior: 'instant'});
                window.scrollBy(0, -margin);
                return true;
            }""",
            [anchor, ANCHOR_MARGIN],
        )
        if not found:
            browser.close()
            # Not SystemExit: main() falls back to a placeholder on Exception,
            # and this is the path a --skip-reports build takes, where the
            # report is itself a placeholder with no such section in it.
            raise RuntimeError(f"no element #{anchor} in {report}")

        page.screenshot(path=str(out))
        browser.close()


def drive_panel(report: Path, out: Path):
    """Record the per-position panel being driven, as an animated WebP.

    A still cannot make this argument. What earns the panel its complexity is
    that a six-position finding inside 398 is a hairline until it is zoomed, and
    that the toolbar will go and find it for you - both of which are motion. So
    the demo is that and nothing else: drag across the flagged region, read one
    position off the tooltip, zoom back out, then let one click of `Next flag`
    land on the finding without the drag.
    """
    import io

    from PIL import Image
    from playwright.sync_api import sync_playwright

    frames, durations = [], []

    with sync_playwright() as playwright:
        browser = playwright.chromium.launch()
        page = browser.new_page(viewport=PANEL_VIEWPORT, device_scale_factor=1)
        page.goto(report.resolve().as_uri())
        page.wait_for_load_state('load')
        page.wait_for_function(CANVAS_PAINTED, timeout=30_000)

        # The forward panel; a report also carries the reversed one, which shows
        # the same finding read from the other end.
        panel = page.locator('.ppv').first
        if not panel.count():
            raise RuntimeError(f"no per-position panel in {report}")
        plot = panel.locator('.ppv-plot')
        panel.scroll_into_view_if_needed()
        page.evaluate(CURSOR_SETUP)

        def shot(ms: int):
            frames.append(Image.open(io.BytesIO(panel.screenshot())).convert('RGB'))
            durations.append(ms)

        def hold(ms: int):
            """Leave the last frame up longer, rather than shooting it twice."""
            durations[-1] = ms

        def move(x: float, y: float, down: bool = False):
            page.mouse.move(x, y)
            page.evaluate('([x, y, d]) => window.__gbqcCursor(x, y, d)',
                          [x, y, down])

        def mapping():
            """Solve the plot's position-to-pixel mapping for the current view.

            Re-solved after every zoom, since it costs two mouse moves.
            """
            box = plot.bounding_box()
            y = box['y'] + box['height'] * 0.35
            xa = box['x'] + box['width'] * 0.35
            xb = box['x'] + box['width'] * 0.75
            page.mouse.move(xa, y)
            pa = panel.evaluate(HOVER_POSITION)
            page.mouse.move(xb, y)
            pb = panel.evaluate(HOVER_POSITION)
            if pa is None or pb is None or pa == pb:
                raise RuntimeError(
                    "could not read the plot's x axis back off the tooltip")
            step = (xb - xa) / (pb - pa)
            return (lambda pos: xa + (pos - pa) * step), y

        def press(action: str, settle_ms: int):
            """Click a toolbar button, with the press visible in the frames."""
            box = panel.locator(f'[data-action="{action}"]').bounding_box()
            cx, cy = box['x'] + box['width'] / 2, box['y'] + box['height'] / 2
            move(cx, cy)
            shot(PANEL_MOVE_MS)
            page.mouse.down()
            move(cx, cy, down=True)
            shot(PANEL_MOVE_MS)
            page.mouse.up()
            move(cx, cy)
            shot(settle_ms)

        flags = panel.evaluate(FLAGGED_POSITIONS)
        if not flags['all']:
            raise RuntimeError(
                f"nothing is flagged in {report}, so there is no finding to "
                f"zoom into - point --report at one that fails a per-position "
                f"check")
        first, last = flags['all'][0], flags['all'][-1]
        centre = (first + last) / 2
        half = max(PANEL_MIN_HALF_WINDOW, (last - first) * PANEL_WINDOW_SPREAD)
        window = (max(1, centre - half), min(flags['end'], centre + half))
        # A Fail if there is one: it is the row whose tooltip has the most in it.
        target = (flags['fails'] or flags['all'])[0]

        # 1. At rest. Hundreds of positions of flat lines, the finding in there
        #    somewhere and invisible - which is the whole reason for the panel.
        shot(PANEL_HOLD_MS)

        # 2. Drag across the flagged region. The selection rectangle grows with
        #    the pointer, then the view snaps to it on release.
        x_of, y = mapping()
        move(x_of(window[0]), y)
        shot(PANEL_MOVE_MS)
        page.mouse.down()
        for i in range(1, PANEL_DRAG_STEPS + 1):
            span = (window[1] - window[0]) * i / PANEL_DRAG_STEPS
            move(x_of(window[0] + span), y, down=True)
            shot(PANEL_MOVE_MS)
        page.mouse.up()
        shot(PANEL_HOLD_MS)

        # 3. Hover onto a flagged position: the tooltip carries the per-class
        #    frequencies and the flag that was raised on them.
        x_of, y = mapping()
        approach = target - (window[1] - window[0]) / 5
        for i in range(1, 6):
            move(x_of(approach + (target - approach) * i / 5), y)
            shot(PANEL_MOVE_MS)
        hold(PANEL_READ_MS)

        # 4. and 5. Back out to the whole sequence, then let the toolbar find
        #    the same place in one click.
        press('reset', PANEL_HOLD_MS)
        press('next', PANEL_END_MS)

        browser.close()

    out.parent.mkdir(parents=True, exist_ok=True)
    # method=6 is the encoder's slowest setting: a few seconds once, against a
    # file every reader of the README downloads.
    frames[0].save(out, save_all=True, append_images=frames[1:],
                   duration=durations, loop=0, quality=72, method=6)
    return frames, durations


def _with_scroll_bar(frame, fraction: float, thumb: int):
    """Draw an overlay scrollbar thumb onto one frame.

    Cheap continuity: it says how far down a long report the frame is, so a
    glide between two sections that happen to look alike reads as movement
    rather than as a cut, and the length of the document is visible from the
    length of the thumb.
    """
    from PIL import Image, ImageDraw

    width, height = frame.size
    right = width - SCROLL_BAR_INSET
    top = round((height - thumb) * fraction)
    overlay = Image.new('RGBA', frame.size, (0, 0, 0, 0))
    ImageDraw.Draw(overlay).rounded_rectangle(
        [right - SCROLL_BAR_WIDTH, top, right, top + thumb],
        radius=SCROLL_BAR_WIDTH // 2, fill=SCROLL_BAR_FILL)
    return Image.alpha_composite(frame.convert('RGBA'), overlay).convert('RGB')


def scroll_tour(report: Path, out: Path):
    """Scroll a whole report top to bottom, pausing on each section.

    What this is for is the scale of the thing: a report is a dozen sections,
    each with a flag, a figure, a plot and a paragraph saying what the figure
    means, and no single frame of it can say so. The stops come from the
    report's own sidebar rather than from a list here, so the tour is however
    long the report is.

    The page is really scrolled rather than cropped out of one tall screenshot,
    which costs a screenshot per frame but is the only way the sidebar behaves:
    it is `position: fixed`, and its scroll-spy moves the highlight down the
    contents as the sections go past.
    """
    import io

    from PIL import Image
    from playwright.sync_api import sync_playwright

    frames, durations = [], []

    with sync_playwright() as playwright:
        browser = playwright.chromium.launch()
        page = browser.new_page(viewport=SCROLL_VIEWPORT, device_scale_factor=1)
        page.goto(report.resolve().as_uri())
        page.wait_for_load_state('load')
        page.add_style_tag(content=PREP_CSS)
        page.wait_for_function(CANVAS_PAINTED, timeout=30_000)

        travel = page.evaluate('Math.max(0, document.documentElement'
                               '.scrollHeight - window.innerHeight)')
        stops = page.evaluate(SCROLL_STOPS, SCROLL_NAV)
        if not stops or not travel:
            raise RuntimeError(
                f"nothing to scroll in {report}: {len(stops)} section(s) in the "
                f"sidebar, {travel}px of travel")

        # (scroll offset, how long that frame is held). A pause is one frame
        # with a long duration, not a dozen identical ones - the same trick
        # --mode panel uses, and what keeps a twelve-second tour to a file this
        # size.
        plan: list[tuple[float, int]] = [(0.0, SCROLL_TOP_MS)]

        def glide_to(target: float, hold_ms: int, step_px: int):
            start = plan[-1][0]
            steps = min(SCROLL_MAX_STEPS, round(abs(target - start) / step_px))
            for step in range(1, steps):
                plan.append((start + (target - start) * step / steps,
                             SCROLL_GLIDE_MS))
            plan.append((target, hold_ms))

        for stop in stops:
            glide_to(min(max(stop['y'] - ANCHOR_MARGIN, 0.0), travel),
                     SCROLL_HOLD_MS, SCROLL_STEP_PX)
        plan[-1] = (plan[-1][0], SCROLL_END_MS)

        # Rewind, then drop the frame that lands back at the top: the first
        # frame already stands there, so the loop closes on it and the wrap is
        # one more step of the rewind rather than a jump cut.
        glide_to(0.0, SCROLL_TOP_MS, SCROLL_RETURN_PX)
        plan.pop()

        # A thumb as tall, against the frame, as one screen is against the whole
        # document - which is what a scrollbar means.
        height = SCROLL_VIEWPORT['height']
        thumb = max(SCROLL_BAR_MIN,
                    round(height * height / (travel + height)))

        for offset, duration in plan:
            page.evaluate('y => window.scrollTo(0, y)', offset)
            page.wait_for_timeout(SCROLL_SETTLE_MS)
            frame = Image.open(io.BytesIO(page.screenshot())).convert('RGB')
            frames.append(_with_scroll_bar(frame, offset / travel, thumb))
            durations.append(duration)

        browser.close()

    print(f"  {len(stops)} section(s), {travel}px of travel")
    out.parent.mkdir(parents=True, exist_ok=True)
    frames[0].save(out, save_all=True, append_images=frames[1:],
                   duration=durations, loop=0, quality=68, method=6)
    return frames, durations


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument('--mode', choices=('still', 'scroll', 'panel'),
                        default='still',
                        help="still: the landing-page shot. scroll: the README's "
                             "tour of a whole report. panel: the README's "
                             "per-position demo.")
    parser.add_argument('--report', type=Path, default=None,
                        help="Report HTML to shoot (default depends on --mode)")
    parser.add_argument('--out', type=Path, default=None,
                        help="Where to write it (default depends on --mode)")
    parser.add_argument('--anchor', default=DEFAULT_ANCHOR,
                        help="Element id to frame the still on")
    parser.add_argument('--require', action='store_true',
                        help="Fail rather than write a placeholder")
    args = parser.parse_args()

    out = args.out or {'still': DEFAULT_OUT, 'scroll': DEFAULT_SCROLL_OUT,
                       'panel': DEFAULT_PANEL_OUT}[args.mode]
    report = args.report or (DEFAULT_SCROLL_REPORT if args.mode == 'scroll'
                             else DEFAULT_REPORT)
    # Nothing in the built site links either animation, so there is no
    # missing-image failure for a placeholder to stand in for and one would
    # reach the README unnoticed. Both fail instead, whatever was asked for.
    require = args.require or args.mode != 'still'

    if not report.exists():
        reason = 'no report to shoot'
        if require:
            raise SystemExit(f"{reason}: {report}")
        print(f"screenshot: {reason} at {report}; writing placeholder")
        write_placeholder(out, reason)
        return 0

    try:
        import playwright.sync_api  # noqa: F401
    except ImportError:
        reason = 'playwright not installed'
        if require:
            raise SystemExit(
                f"{reason}. Install it with `pip install playwright && "
                f"python -m playwright install chromium`, or drop --require "
                f"to accept a placeholder."
            ) from None
        print(f"screenshot: {reason}; writing placeholder")
        write_placeholder(out, reason)
        return 0

    try:
        if args.mode == 'scroll':
            frames, durations = scroll_tour(report, out)
        elif args.mode == 'panel':
            frames, durations = drive_panel(report, out)
        else:
            capture(report, out, args.anchor)
    except Exception as error:
        # Most often the Chromium build is missing while the Python package is
        # present, which is its own error class rather than an ImportError.
        if require:
            raise
        print(f"screenshot: {type(error).__name__}: {error}; "
              f"writing placeholder")
        write_placeholder(out, 'browser unavailable')
        return 0

    size = out.stat().st_size
    if args.mode == 'still':
        print(f"screenshot: {out.relative_to(ROOT)} "
              f"({size / 1024:.0f} KiB) from {report.relative_to(ROOT)}"
              f"#{args.anchor}")
    else:
        print(f"{args.mode} animation: {out.relative_to(ROOT)} "
              f"({size / 1024:.0f} KiB, {len(frames)} frames, "
              f"{sum(durations) / 1000:.1f}s) from {report.relative_to(ROOT)}")
    return 0


if __name__ == '__main__':
    sys.exit(main())

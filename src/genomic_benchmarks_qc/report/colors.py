"""The flag colors, in one place.

The plots draw them with matplotlib and the report draws them with CSS, so they
cannot live in either one: a divergence would mean a position shaded amber in a
figure and red in the sidebar beside it.
"""

FAIL_COLOR = '#c62828'   # .status-fail background-color
WARN_COLOR = '#f57f17'   # .status-warn background-color
PASS_COLOR = '#2e7d32'   # .status-pass background-color
# Unknown is not a fourth severity - it says the check was not run - so it is
# grey rather than another point on the red-to-green scale.
UNKNOWN_COLOR = '#8a8983'  # .status-unknown background-color

# The two classes, in every plot and in the interactive figure's legend. Kept
# here rather than in the plotting module so the canvas and matplotlib draw the
# same two blues.
CLASS_COLORS = ('#003D99', '#66A3FF')

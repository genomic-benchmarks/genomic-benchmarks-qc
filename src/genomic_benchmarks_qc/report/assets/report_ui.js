/* Behaviour shared by the report pages: the explanation toggles, the table of
 * duplicate sequences, the nav's current-section highlight, and the
 * needs-attention filter.
 *
 * Written as plain ES5 with no dependencies, because a report is a file on
 * someone's disk that may be opened years from now, in whatever browser is to
 * hand, possibly from a share over file://.
 */
(function () {
  'use strict';

  /* ---------------- explanations ----------------
   * Called from onclick attributes in the markup, so it has to be global. */
  window.toggleExplanation = function (elementId) {
    var element = document.getElementById(elementId);
    if (element) element.classList.toggle('visible');
  };

  /* ---------------- duplicate sequences ----------------
   * The list arrives as JSON in a script element rather than as a JS literal,
   * so a sequence containing markup cannot break out of the script. */
  function fillDuplicateSequences() {
    var body = document.querySelector('#sequence-duplication-levels tbody');
    var data = document.getElementById('duplicate-sequences');
    if (!body || !data) return;
    var sequences;
    try { sequences = JSON.parse(data.textContent); } catch (e) { return; }
    sequences.forEach(function (sequence) {
      var row = document.createElement('tr');
      var cell = document.createElement('td');
      cell.textContent = sequence;
      cell.className = 'sequence_column';
      row.appendChild(cell);
      body.appendChild(row);
    });
  }

  /* ---------------- state ----------------
   * Copy the verdict each check already shows onto the section and its nav
   * entry, so the filter and the stylesheet can select on it. */
  var STATES = /status-(pass|warn|fail|unknown)/;

  function stateOf(el) {
    var icon = el.querySelector('.status-icon');
    var match = icon && icon.className.match(STATES);
    return match ? match[1] : null;
  }

  function stampStates() {
    var byId = {};
    var counts = { pass: 0, warn: 0, fail: 0, unknown: 0 };
    [].slice.call(document.querySelectorAll('.content section')).forEach(function (section) {
      var state = stateOf(section);
      if (!state) return;
      section.setAttribute('data-state', state);
      counts[state]++;
      if (section.id) byId[section.id] = state;
    });
    [].slice.call(document.querySelectorAll('.sidebar > .sidebar-item')).forEach(function (item) {
      var link = item.querySelector('a[href^="#"]');
      var state = link && byId[link.getAttribute('href').slice(1)];
      if (state) item.setAttribute('data-state', state);
    });
    return counts;
  }

  /* ---------------- which check am I looking at ----------------
   * Exactly one nav entry is marked. An IntersectionObserver leaves none marked
   * in the gaps between sections and two marked where they overlap, which makes
   * the nav look broken. */
  function markCurrentSection() {
    var rows = [];
    [].slice.call(document.querySelectorAll('.sidebar > .sidebar-item')).forEach(function (item) {
      var link = item.querySelector('a[href^="#"]');
      var target = link && document.getElementById(link.getAttribute('href').slice(1));
      if (target) rows.push({ item: item, target: target });
    });
    if (!rows.length) return;

    var queued = false;
    function mark() {
      queued = false;
      // The current check is the last one whose top has passed the top of the
      // viewport. A line further down the page marks whatever sits a third of
      // the way in, which reads as the wrong entry when the check you are
      // looking at is the one at the top; the few pixels of slack only absorb
      // sub-pixel scroll positions after clicking a nav link.
      var line = 8;
      var active = rows[0];
      rows.forEach(function (row) {
        if (row.target.getBoundingClientRect().top <= line) active = row;
      });
      // at the very bottom the last section wins, wherever it starts
      if (window.innerHeight + window.pageYOffset >= document.body.scrollHeight - 4) {
        active = rows[rows.length - 1];
      }
      rows.forEach(function (row) {
        row.item.classList.toggle('is-active', row === active);
      });
    }
    function schedule() {
      if (queued) return;
      queued = true;
      if (window.requestAnimationFrame) window.requestAnimationFrame(mark);
      else setTimeout(mark, 60);
    }
    window.addEventListener('scroll', schedule, { passive: true });
    window.addEventListener('resize', schedule);
    mark();
  }

  /* ---------------- the needs-attention filter ---------------- */
  var MIN_CHECKS_FOR_FILTER = 3;   // below that, the whole report fits on a screen

  function addFilter(counts) {
    var flagged = counts.fail + counts.warn;
    var total = flagged + counts.pass + counts.unknown;
    if (total < MIN_CHECKS_FOR_FILTER) return;

    var box = document.createElement('div');
    box.className = 'qc-filter';
    box.setAttribute('role', 'group');
    box.setAttribute('aria-label', 'Which checks to show');
    box.innerHTML =
      '<button type="button" data-mode="attention" aria-pressed="false">Needs attention'
      + (flagged ? ' (' + flagged + ')' : '') + '</button>'
      + '<button type="button" data-mode="all" aria-pressed="true">All (' + total + ')</button>';

    var buttons = [].slice.call(box.querySelectorAll('button'));
    buttons.forEach(function (button) {
      button.addEventListener('click', function () {
        var only = button.getAttribute('data-mode') === 'attention';
        document.body.classList.toggle('qc-attention', only);
        buttons.forEach(function (other) {
          other.setAttribute('aria-pressed', String(other === button));
        });
        // a hidden .ppv-plot measures 0 wide, so the figure re-measures on show
        window.dispatchEvent(new Event('resize'));
      });
    });

    var sidebar = document.querySelector('.sidebar');
    if (!sidebar) return;
    sidebar.insertBefore(box, sidebar.querySelector('h2') || sidebar.firstChild);
  }

  fillDuplicateSequences();
  addFilter(stampStates());
  markCurrentSection();
}());

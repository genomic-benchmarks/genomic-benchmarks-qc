/* Behaviour shared by the report pages: the explanation toggles and the table
 * of duplicate sequences.
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

  fillDuplicateSequences();
}());

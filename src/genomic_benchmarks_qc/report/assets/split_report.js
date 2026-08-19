/* Behaviour specific to the split report: showing one alignment at a time.
 * Called from an onclick attribute in the markup, so it has to be global. */
window.toggleAlignment = function (id, button) {
  var row = document.getElementById(id);
  if (!row) return;
  var hidden = row.style.display === 'none';
  row.style.display = hidden ? 'table-row' : 'none';
  button.textContent = hidden ? 'Hide' : 'Show';
};

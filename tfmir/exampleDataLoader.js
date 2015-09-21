function loadExampleData() {
	console.log('Mimic sample data upload');
	$('#mirnaDemo').attr('value', true);
	$('#mrnaDemo').attr('value', true);

	uploadFile('miRNA');
	uploadFile('mRNA');	

	$('#mirnaDemo').attr('value', false);
	$('#mrnaDemo').attr('value', false);

	$('select[id=disease] option[selected="selected"]').removeAttr('selected');
	
	//$('select[id=disease] option[value="Breast Neoplasms"]').attr('selected','selected');
	
	var text1 = 'Breast Neoplasms';


	$("select[id=disease] option").filter(function()
	{
		return $(this).hasClass('highlight');
	}).removeClass('highlight').removeAttr('selected');

	$("select[id=disease] option").filter(function() {
    //may want to use $.trim in here
    	return $(this).text() == text1; 
	}).attr('selected', true).addClass('highlight');

	//$('#disease').chosen().change();
	$("#disease").trigger("chosen:updated");

	$("select[id=disease] option").filter(function()
	{
		return $(this).hasClass('highlight');
	}).removeClass('highlight').removeAttr('selected');

	$("select[id=evidence] option").filter(function() {
    //may want to use $.trim in here
    	return $(this).text() == 'Both'; 
	}).attr('selected', true).addClass('highlight');

	//$('#disease').chosen().change();
	$("#evidence").trigger("chosen:updated");
	$('#disease').val("Breast Neoplasms");
	$('#evidence').val("Both");

	log("Example files for breast cancer have been loaded into your session. Check your p-Value and experimental evidence and click the processing button to start analysis.")
}
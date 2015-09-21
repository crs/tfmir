<?php header('Content-Type: text/javascript; charset=UTF-8'); ?>

$('#motifs').change(function() {
	console.log(this.text());
});

if (typeof String.prototype.startsWith != 'function') {
  // see below for better implementation!
  String.prototype.startsWith = function (str){
    return this.indexOf(str) == 0;
  };
}

function highlightMotif(value) {
  var arr = value.split(",");
  console.log(arr);

  var cy = $('.cytoPanel').cytoscape('get');

  cy.elements().unselect();

  console.log('[id="' + arr[1] + '"]');
  cy.nodes('[id="' + arr[0] + '"]').select();
  cy.nodes('[id="' + arr[1] + '"]').select();
  cy.nodes('[id="' + arr[2] + '"]').select();
}

function extendSelectedNeighbours() {
  var cy = $('.cytoPanel').cytoscape('get');

  if (cy.nodes(':selected').length < 1) {
    alert("No nodes selected. Please select a node to start with.");
  }

  cy.nodes(':selected').neighborhood().select();
}

function addSelector(motifSearchResponse) {
	$('#motifs').css('visibility','visible');
	//$('#motifs').innerHTML = '';
//  var cy = $(".cytoPanel").cytoscape('get');

	$.each(motifSearchResponse.Motifs, function(i, val) {
  		//$("#" + i).append(document.createTextNode(" - " + val));
  		//$('#motifs').append('<option onclick="$(\'.cytoPanel\').cytoscape(\'get\').nodes(\'[id="'+val.tf+'"]\'); console.log(\'blubb\'); ">' + val.type + ": " + val.tf + "-" + val.mirna + " " + val.gene + "</option>");
      $('#motifs').append('<option value="' + val.tf + ',' + val.mirna + ',' + val.gene + '">' + val.type + ": " + val.tf + "-" + val.mirna + " " + val.gene + "</option>");
      //$('#motifs').append('<option>' + val.gene + "</option>");
	});
}

function motifSearch(file) {
	var xmlhttp;
	console.log("Start motif search");
  $('#motifloader').css('visibility','visible')

	if (window.XMLHttpRequest) {// code for IE7+, Firefox, Chrome, Opera, Safari
  		xmlhttp=new XMLHttpRequest();
 	} else {// code for IE6, IE5
  		xmlhttp=new ActiveXObject("Microsoft.XMLHTTP");
  }
	
	xmlhttp.onreadystatechange=function()
  	{
      console.log(xmlhttp);
  		if (xmlhttp.readyState==4 && xmlhttp.status==200) {
    		//document.getElementById("myDiv").innerHTML=xmlhttp.responseText;

        if (xmlhttp.responseText.startsWith("<")) {
          alert("There has been an error sending the request. Probably cytoscape is down. :( Please let us know.");
        } else {
          var motifSearchResponse = JSON.parse(xmlhttp.responseText);
        //console.log(xmlhttp);
        //console.log(xmlhttp.responseText);

          addSelector(motifSearchResponse);
        }
        $('#motifloader').css('visibility','hidden');
        $('#startMotifSearchButton').remove();
    	} //else if (xmlhttp.readyState==0) {
        //alert("Could not perform motif search. Maybe cytoscape is offline?");
        //alert(xmlhttp.status);
      //}
  	}

	//xmlhttp.open("GET","motifSearch.php?dataset="+file, true);
  xmlhttp.open("GET","motifSearchPHP.php?dataset="+file, true);
	xmlhttp.send();
}
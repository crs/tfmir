<?php header('Content-Type: text/javascript; charset=UTF-8'); ?>


function showCos(folder, dataset, type) {
	var tf_mirna_gene = ($("#motifs").val()).split(",");

	var tf = tf_mirna_gene[0];
	var mirna = tf_mirna_gene[1];
	var gene = tf_mirna_gene[2];

    //cy.elements('node#'+tf);

    window.open('motif.php?folder='+folder+'&dataset='+dataset+'&tf='+tf+'&mirna='+mirna+'&gene='+gene+'&type='+type);

  //var generator=window.open('','Motif view');

  //generator.document.write('<html><head><title>Popup</title></head><body>');
  //generator.document.write('<div id=motifNet></div>');
  //generator.document.write('<script type="text/javascript">var json = ' + cy.json() +'</script>');
  //generator.document.write('<pre>'+cy.json()+'</pre>');
  //generator.document.write('<script type="text/javascript" src="motifNet.js.php"></script>');
  //generator.document.write('</body></html>');


  //generator.document.close();

}

function calculateFunctionalSimilarity(folder,genes,motif,co) {
	//$.get('functionalSimilarity.php?folder='+folder+'&genes='+genes+'&motif='+motif, function(data) {
	//	console.log(data);
	//});
	var xmlhttp;
	console.log("Start functional similarity calculation");
	
	console.log($('#similarityLoader'));
	$('#similarityLoader').css('visibility','visible');

	if (window.XMLHttpRequest) {// code for IE7+, Firefox, Chrome, Opera, Safari
		xmlhttp=new XMLHttpRequest();
 	} else {// code for IE6, IE5
 		xmlhttp=new ActiveXObject("Microsoft.XMLHTTP");
 	}

 	xmlhttp.onreadystatechange=function() {
 		console.log(xmlhttp);
 		
 		if (xmlhttp.readyState==4 && xmlhttp.status==200) {
    		//document.getElementById("myDiv").innerHTML=xmlhttp.responseText;

    		if (xmlhttp.responseText.startsWith("0")) {
    			alert("There has been an error sending the request. Probably R is not connectable. :( Please let us know.");
    			$('#similarityLoader').css('visibility','hidden');
    		} else {
    			//var motifSearchResponse = console.log(xmlhttp.responseText);

    			var img = new Image();

    			$(img).load(function() {
    				$("#funcSimImage").append(this);
    				$(this).fadeIn();
    			}).error(function() { alert("Image has not been created. Sorry there has been an error.") })
    			.attr('src', "uploads/" + folder+"/"+co+'-'+motif+'.png');
    			$('#similarityLoader').css('visibility','hidden');
    			$('#calculationButtonFS').remove();
    		}

    	} //else { alert("Sorry, server was not available"); }
    	//window.location.reload();
    	

    } //else if (xmlhttp.readyState==0) {
        //alert("Could not perform motif search. Maybe cytoscape is offline?");
        //alert(xmlhttp.status);
      //}
  console.log(genes);
  console.log(motif);
  console.log(co);
  xmlhttp.open("GET",'functionalSimilarity.php?folder='+folder+'&genes='+genes+'&motif='+motif+'&type='+co, true);
  xmlhttp.send();
}
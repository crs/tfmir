<?php header('Content-Type: text/javascript; charset=UTF-8'); 
session_start();
 $url = '//'.$_SERVER['HTTP_HOST'].dirname($_SERVER['PHP_SELF']);

?>


function showDominatingSet(file) {
	$('#domsetloader').css('visibility','hidden');
	$("#dialog").dialog({
    			autoOpen: false,
			    position: 'center' ,
			    title: 'Dominating Nodes',
			    draggable: false,
			    resizable : false,
			    modal : true,
			});

	var domsetFile = './uploads/<?php echo session_id() ?>/' + file + '/dominatingSet.txt'
	var loc  = '<?php echo $url; ?>/exists.php?fileName='+domsetFile;
	console.log(loc);
	$.get(loc, function(data) {
		if (data == 1) {
			console.log(data);
			//$("#dialog code").load(domsetFile, function() {
				//$("#dialog").dialog("open");
			//	console.log(domsetFile);
			//});

			$.get(domsetFile, function(data) {
				var str = data.split("\n");
				console.log(str);
				var cy = $('.cytoPanel').cytoscape('get');

				cy.elements().unselect();
				for (var i = 0; i < str.length; i++) {

					if (str[i].length == 0) continue;
					console.log(str[i]);
					cy.nodes('[id="' + str[i] + '"]').select();
					cy.nodes('[id="' + str[i].toLowerCase() + '"]').select();
				}
			});
		} else {
			$('#domsetloader').css('visibility','visible')
			$.get('<?php echo $url; ?>/cgi-bin/domSet.cgi?session=<?php echo session_id(); ?>&interactions=' + file, function(domsetcompute) {
				console.log("Data not found, try to compute");
				console.log(domsetcompute);
				if (domsetcompute == 1) {
					console.log("successfully computed dominating set");
					showDominatingSet(file);
				} else {
					alert('Could not compute dominating set');
				}
			});
		}
	});
}
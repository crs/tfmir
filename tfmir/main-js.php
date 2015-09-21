<?php header('Content-Type: text/javascript; charset=UTF-8'); 
//$url = (!empty($_SERVER['HTTPS'])) ? "https://".$_SERVER['SERVER_NAME'].$_SERVER['REQUEST_URI'] : "http://".$_SERVER['SERVER_NAME'].$_SERVER['REQUEST_URI']; 
 $url = '//'.$_SERVER['HTTP_HOST'].dirname($_SERVER['PHP_SELF']);
 ?>
//if (!serverPath) {
//   var serverPath = window.location.protocol+"//"+window.location.host+"/tfmir/"+;
//   console.log("Server Path: " + serverPath);
//}

//var url = "<?php echo $url ?>";
var scriptEls = document.getElementsByTagName( 'script' );
var thisScriptEl = scriptEls[scriptEls.length - 1];
var scriptPath = thisScriptEl.src;
var scriptFolder = scriptPath.substr(0, scriptPath.lastIndexOf( '/' )+1 );

//console.log( [scriptPath, scriptFolder] );

var serverPath = '<?php echo $url; ?>';

function log(message) {
	document.getElementById('entries').innerHTML = "[" + formatDate(new Date(), "HH:mm:ss") + "] " + message + "<br/>"+document.getElementById('entries').innerHTML;
	console.log("[" + formatDate(new Date(), "HH:mm:ss") + "] " + message);
	//document.getElementById('entries').innerHTML =
}

function baseName(str) {
   return str.split(/[\\/]/).pop();
}

function setCookie(event) {
	console.log([event.type, event.value]);
	console.log([event]);
	if (event.type == 'checkbox') {
			console.log("Setting cookie for " + event.id+'='+event.checked);
			document.cookie = event.target.name+'='+event.target.checked;
	}
	else {
		console.log(event.value);
		console.log("Setting cookie for " + event.id+'='+event.value);
		document.cookie = event.id+'='+event.value ;
	}

}

function checkRange(lower, upper) {
	var val = parseFloat(event.target.value);
	if ((val < lower) || (val > upper)) {
		log("Value \"" + val + "\" not allowed, has to be in range between "+lower+" and "+upper);
		console.log("Error: value not in allowed range");
		//event.target.style.background.color="#FF0000";
		event.target.value=readCookie(event.target.name);
		return false;
	}
	return true;
}

function readCookie(name) {
	var cookiename = name + "=";
	var ca = document.cookie.split(';');
	for(var i=0;i < ca.length;i++)
	{
		var c = ca[i];
		while (c.charAt(0)==' ') c = c.substring(1,c.length);
		if (c.indexOf(cookiename) === 0) return c.substring(cookiename.length,c.length);
	}
	return null;
}

function fileSelected(event) {
//alert(button);
	console.log(event);
	for (var i = 0; i < event.files.length; i++) { console.log(event.files[i]); };
  var file = event.files[0];//document.getElementById('fileToUpload').files[0];
  var button = event.id + '-button';
  console.log('Button id:' + button);
  log(file);
	if (file) {
		if ((Math.round(file.size * 100 / (1024 * 1024)) / 100) > 128) {
			log(file.name + " is too large, sorry.");
			//document.getElementById(button).setAttribute('disabled','disabled');
		}	else {
			var fileSize = 0;
			if (file.size > 1024 * 1024)
				fileSize = (Math.round(file.size * 100 / (1024 * 1024)) / 100).toString() + 'MB';
			else
				fileSize = (Math.round(file.size * 100 / 1024) / 100).toString() + 'KB';

			log(event.id + " file: " + file.name + " " + fileSize + " " + file.type);
			//document.getElementById(button).removeAttribute('disabled');
    //document.getElementById('fileName').innerHTML = 'Name: ' + file.name;
    //document.getElementById('fileSize').innerHTML = 'Size: ' + fileSize;
    //document.getElementById('fileType').innerHTML = 'Type: ' + file.type;
			}
		} else {
			log('No file could be loaded');
		}
	}

function getmiRNAResultTable() {
	var xhr = new XMLHttpRequest();

	xhr.addEventListener("load", miRNAResultTableLoaded, false);

	var parameter = 'miRNATableFile=./uploads/' + document.getElementById('session').value + '/transmir.res.txt';
	var url = serverPath + 'tableFunctions.php';
	console.log(url);
	console.log(parameter);

	xhr.open('POST', url);
	xhr.setRequestHeader("Content-type", "application/x-www-form-urlencoded");
	//xhr.setRequestHeader("Content-length", parameter.length);
	//xhr.setRequestHeader("Connection", "close");

	xhr.send(parameter);
}

function miRNAResultTableLoaded(evt) {
	document.getElementById('resultBox').innerHTML = evt.target.response;
	$('#miRNATable').DataTable();
}


function checkResults() {
	//var form = document.getElementById('execute');
	//var fd = new FormData(form);

	var sFolderPath = "uploads/" + document.getElementById('session').value;
	//var fso = new ActiveXObject("Scripting.FileSystemObject");

	//if (!fso.FolderExists(sFolderPath)) {
	//	alert('Folder does not exist!');
	//}

	var resultFolders = ['tf-gene', 'tf-mirna', 'mirna-mirna', 'mirna-gene', 'all', 'disease',  'mRNA'];

	for (i = 0; i < resultFolders.length; i++) {
			//alert(resultFolders[i]);
			setActiveIfAvailable(sFolderPath, resultFolders[i]);
	}
}

// get output from uri, should not take to long as this is an synchronous call
function setActiveIfAvailable(url, element) {
	var urlData;
	var loc;
	if (element == 'mRNA') {
		loc = '<?php echo $url; ?>/exists.php?fileName=' + url + '/' + element;
	} else {
		loc = '<?php echo $url; ?>/exists.php?fileName=' + url + '/' + element + '/summary.txt';
	}
	console.log(loc);
	$.get(loc, function(data) {
									var selector = '.resultButton[id=\'' + element + '\']';
									var selector = 'div #' + element + ':not(select)';

									if (element == 'disease') {
										var selector = 'div .diseaseButton';
									}

									console.log($(selector));	

									if (data == 1) {
									if (element == 'mRNA') {
										selector = '#processingButton';
										if ($(selector).hasClass('inactive')) $(selector).removeClass('inactive');
										document.getElementById('processingButton').onclick = function () {
  											startProcessing(); checkResults(); return false;};
  											$('#processingButton').attr('title',"Ready! Click to start processing.");
										return;
									}

									$('.hiddenResultPanel').css('display','block');
									var session = document.getElementById('session').value;
									
									if ($(selector).hasClass('inactive')) {
										$(selector).removeClass('inactive'); 
									}
										var elementPath = url + '/' + element;
									//$(selector).css('background-image', 'url('+url+'/'+element+'/venn.png)');
										
										var imagePath;
										if (element == 'all' | element == 'disease') {
											//imagePath = url + '/' + element + '/degree.png';
										} else {
											imagePath = url + '/' + element + '/venn.png';
										}

										if (imagePath) {
										$(selector).css({
											'background-image':'url(' + imagePath + ')',
											'background-position':'50%',
											'background-size':'contain',
											'background-repeat':'no-repeat'
										});
										}
										
										

										// onclick attempt, did not work
										//$(selector).onclick = "function() { window.open('www.google.de')});"

										$(selector).html('<a href="viewResult.php?folder=' + session + '&dataset='+ element +'">' + $(selector).html() + '</a>');
										
										var text = $.get('resultButtonSummary.php?folder=' + session + '&dataset='+ element, function(data) { $(selector).append('<div class="buttonInfo">'+data+'</div>');});
									} else {
									//$(selector).addClass('inactive');
								} 
							});
}

function startProcessing() {
	var xhr = new XMLHttpRequest();

	var form = document.getElementById('execute');
	var fd = new FormData(form);

	xhr.upload.addEventListener("progress", processingProgress, false);
	xhr.addEventListener("load", processingComplete, false);
  xhr.addEventListener("error", uploadFailed, false);
	xhr.addEventListener("abort", uploadCanceled, false);


	
	//$("#loading_dialog").loading();
	//$("#processingButton").loading();
	$('#pageframe').css('background-color','rgba(0,0,0,0.8)');
	$('#processingButton').css({'background-color':'#fffb7d', 'background-image':'url(img/loader.gif)'});
	$('#processingButton').html('<h4>Please stand by...</h4>');
	$('#processingButton').attr('title','Please wait until the processing is complete.');
	//xhr.open('POST', "<?php echo $url; ?>/cgi-bin/execute.cgi", true);
	xhr.open('POST', "<?php echo $url; ?>/executeTFmir.php", true);
	xhr.send(fd);

	//$("#loading_dialog").loading("loadStart","Please stand by...");
}

function processingComplete(evt) {
	//log("Processing complete! <a href='' onclick='getmiRNAResultTable(); return false;'>Click to show</a>");
	//$('#results').css('display','block');
	log("Job finished");
	//$("#processingButton").loading("loadStop");
	$('#processingButton').css({'background-color':'#66d15c', 'background-image':'none'});
	$('#processingButton').html('<h4>Processing Complete!</br><small>Explore the networks by clicking the buttons below or click to restart</small></h4>');
	$('#processingButton').attr('title','The are results ready. Click to start over the processing job!');
	checkResults();
}

function exportImage() {
	window.open($('.cytoPanel').cytoscape('get').png(), 'Network View');
}

function processingProgress(evt) {
	log("Job now running!");
}

function uploadFile(formId) {
	document.getElementById('meter').style.visibility='visible';

  var xhr = new XMLHttpRequest();
  
  var form = document.getElementById(formId);
  
  var fd;
  if (document.getElementById('mirnaDemo').value == 'true') {
  	console.log("Loading Demo");
  	fd = createFormData(formId);
  	xhr.filename = fd.filename;
  } else {
  	console.log("Uploading files from form");
  	fd = new FormData(form);
  	xhr.filename = form.uploadedfile.value;
  }
  

  /* event listners */
  xhr.upload.addEventListener("progress", uploadProgress, false);
  xhr.addEventListener("load", uploadComplete, false);
  xhr.addEventListener("load", function() {
  	if (formId === 'mRNA') {
  		$('#processingButton').removeClass('inactive');
  		document.getElementById('processingButton').onclick = function () {
  			startProcessing(); checkResults(); return false;};
  		$('#processingButton').attr('title',"Ready! Click to start processing.");
  		log('mRNA received, ready to start processing');
  	}	
  });
  
  xhr.addEventListener("error", uploadFailed, false);
  xhr.addEventListener("abort", uploadCanceled, false);
  xhr.formtype = form.type.value;
	
  /* Be sure to change the url below to the url of your upload server side script */
  xhr.open("POST", "<?php echo $url; ?>/upload.php");
  xhr.send(fd);
	document.getElementById(formId + '-filename').value = xhr.filename;
	
	console.log("FormData filename: " + fd.filename + "ID: " + "#current_" + formId);
	//log(document.getElementById(formId + '-filename').value);
	if (formId === 'mRNAfileToUpload') {
		log('mRNA')
	}
}

function uploadProgress(evt) {
  if (evt.lengthComputable) {
    var percentComplete = Math.round(evt.loaded * 100 / evt.total);
    //document.getElementById('progressNumber').innerHTML = percentComplete.toString() + '%';
		document.getElementById('progressSpan').style.width = percentComplete + '%';
  }
  else {
    document.getElementById('progressNumber').innerHTML = 'unable to compute';
  }
}

function uploadComplete(evt) {
  /* This event is raised when the server send back a response */
  //document.getElementById('progressNumber').innerHTML = evt.target.responseText;
	document.getElementById('meter').style.visibility='hidden';
	log(evt.target.response);
	document.getElementById('current_' + this.formtype).innerHTML = baseName(this.filename);
}

function uploadFailed(evt) {
  alert("There was an error attempting to upload the file.");
	document.getElementById('meter').style.visibility='hidden';
}

function uploadCanceled(evt) {
  alert("The upload has been canceled by the user or the browser dropped the connection.");
	log(evt.target.response);
	document.getElementById('meter').style.visibility='hidden';
	console.log(this.status);
}


(function($) {
    $.widget("artistan.loading", $.ui.dialog, {
        options: {
            // your options
            spinnerClassSuffix: 'spinner',
            spinnerHtml: 'Loading',// allow for spans with callback for timeout...
            maxHeight: false,
            maxWidth: false,
            minHeight: 80,
            minWidth: 220,
            height: 80,
            width: 220,
            modal: true,
        },

        _create: function() {
            $.ui.dialog.prototype._create.apply(this);
            // constructor
            $(this.uiDialog).children('*').hide();
            var self = this,
            options = self.options;
            self.uiDialogSpinner = $('.ui-dialog-content',self.uiDialog)
                .html(options.spinnerHtml)
                .addClass('ui-dialog-'+options.spinnerClassSuffix);
        },
        _setOption: function(key, value) {
            var original = value;
            $.ui.dialog.prototype._setOption.apply(this, arguments);
            // process the setting of options
            var self = this;

            switch (key) {
                case "innerHeight":
                    // remove old class and add the new one.
                    self.uiDialogSpinner.height(value);
                    break;
                case "spinnerClassSuffix":
                    // remove old class and add the new one.
                    self.uiDialogSpinner.removeClass('ui-dialog-'+original).addClass('ui-dialog-'+value);
                    break;
                case "spinnerHtml":
                    // convert whatever was passed in to a string, for html() to not throw up
                    self.uiDialogSpinner.html("" + (value || '&#160;'));
                    break;
            }
        },
        _size: function() {
            $.ui.dialog.prototype._size.apply(this, arguments);
        },
        // other methods
        loadStart: function(newHtml){
            if(typeof(newHtml)!='undefined'){
                this._setOption('spinnerHtml',newHtml);
            }
            this.open();
        },
        loadStop: function(){
            this._setOption('spinnerHtml',this.options.spinnerHtml);
            this.close();
        }
    });
})(jQuery);


function createFormData(formId) {

	var formData = new FormData();

	formData.type = 'POST';
	formData.append("type", formId);
	var content = getExampleFile(formId);
	var blob = new Blob([content], { type: "text/plain"});

	formData.append("uploadedfile", blob, formId + '.bc.txt');
	formData.filename = formId + '.bc.txt';
	console.log(formData.filename);
	console.log("Created FormData");
	return(formData);
}

function getExampleFile(formId) {
	if (formId == "miRNA") {
		return <?php echo json_encode(file_get_contents("./mirna.bc.txt")); ?>;
	} else if (formId == "mRNA") {
		return <?php echo json_encode(file_get_contents("./mrna.bc.txt")); ?>;
	}
}

function msieversion() {

        var ua = window.navigator.userAgent;
        var msie = ua.indexOf("MSIE ");

        return (msie > 0 || !!navigator.userAgent.match(/Trident.*rv\:11\./));      // If Internet Explorer, return version number
            //alert(parseInt(ua.substring(msie + 5, ua.indexOf(".", msie))));
}

function clearSession() {
	if (confirm('This will remove all generated results. Are you sure?')) {
		console.log('delete');
		$.get('deleteSession.php', function(data) {
			log('Deleted');
			location.reload();
		});
	} else {
		console.log('keep');
	}
}

//if (!serverPath) {
//   var serverPath = window.location.protocol+"//"+window.location.host+"/tfmir/"+;
//   console.log("Server Path: " + serverPath);
//}

var scriptEls = document.getElementsByTagName( 'script' );
var thisScriptEl = scriptEls[scriptEls.length - 1];
var scriptPath = thisScriptEl.src;
var scriptFolder = scriptPath.substr(0, scriptPath.lastIndexOf( '/' )+1 );

console.log( [scriptPath, scriptFolder] );

var serverPath = scriptFolder;

function log(message) {
	document.getElementById('entries').innerHTML = "[" + formatDate(new Date(), "HH:mm:ss") + "] " + message + "<br/>"+document.getElementById('entries').innerHTML;
	console.log("[" + formatDate(new Date(), "HH:mm:ss") + "] " + message);
	//document.getElementById('entries').innerHTML =
	// 'FACKTHIS';
}

function baseName(str) {
   return str.split(/[\\/]/).pop();
}

function setCookie(event) {
	if (event.target.type == 'checkbox') {
			console.log("Setting cookie for " + event.target.name+'='+event.target.checked);
			document.cookie = event.target.name+'='+event.target.checked;
	}
	else {
		console.log(event.target.value);
		console.log("Setting cookie for " + event.target.name+'='+event.target.value);
		document.cookie = event.target.name+'='+event.target.value ;
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

	var resultFolders = ['tf-gene', 'tf-mirna', 'mirna-mirna', 'mirna-gene', 'all', 'mRNA'];

	for (i = 0; i < resultFolders.length; i++) {
			//alert(resultFolders[i]);
			setActiveIfAvailable(sFolderPath, resultFolders[i]);
	}

}

// get output from uri, should not take to long as this is an synchronous call
function setActiveIfAvailable(url, element) {
	var urlData;
	var loc = window.location.protocol + '//' + window.location.hostname + '/tf/exists.php?fileName=' + url + '/' + element;
	console.log(loc);
	$.get(loc, function(data) {
									if (data == 1) {
									var selector = '#' + element;

									if (element == 'mRNA') {
										selector = '#processingButton';
										if ($(selector).hasClass('inactive')) $(selector).removeClass('inactive');
										return;
									}
									var session = document.getElementById('session').value;
									
									if ($(selector).hasClass('inactive')) {
										$(selector).removeClass('inactive'); 
									
										var elementPath = url + '/' + element;
									//$(selector).css('background-image', 'url('+url+'/'+element+'/venn.png)');
										
										var imagePath;
										if (element == 'all') {
											imagePath = url + '/' + element + '/degree.png';
										} else {
											imagePath = url + '/' + element + '/venn.png';
										}
										$(selector).css({
											'background-image':'url(' + imagePath + ')',
											'background-position':'center',
											'background-size':'contain',
											'background-repeat':'no-repeat'
										});
										
										

										// onclick attempt, did not work
										//$(selector).onclick = "function() { window.open('www.google.de')});"

										$(selector).html('<a href="viewResult.php?folder=' + session + '&dataset='+ element +'">' + $(selector).html() + '</a>');
										
										var text = $.get('resultButtonSummary.php?folder=' + session + '&dataset='+ element, function(data) { $(selector).append('<span>'+data+'</span>');});
									}
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

	xhr.open('POST', "http://localhost/tf/cgi-bin/execute.cgi", true);
	xhr.send(fd);
}

function processingComplete(evt) {
	//log("Processing complete! <a href='' onclick='getmiRNAResultTable(); return false;'>Click to show</a>");
	//$('#results').css('display','block');
	//alert(evt.target.response);
	log("Job finished");
	checkResults();
}

function processingProgress(evt) {
	log("Job now running!");
}

function uploadFile(formId) {
	document.getElementById('meter').style.visibility='visible';

  var xhr = new XMLHttpRequest();
  var form = document.getElementById(formId);
	var fd = new FormData(form);

  /* event listners */
  xhr.upload.addEventListener("progress", uploadProgress, false);
  xhr.addEventListener("load", uploadComplete, false);
  xhr.addEventListener("load", function() {
  	if (formId === 'mRNA') {
  		$('#processingButton').removeClass('inactive');
  		log('mRNA received, ready to start processing');
  	}	
  });
  
  xhr.addEventListener("error", uploadFailed, false);
  xhr.addEventListener("abort", uploadCanceled, false);
	xhr.formtype = form.type.value;
	xhr.filename = form.uploadedfile.value;
  /* Be sure to change the url below to the url of your upload server side script */
  xhr.open("POST", "http://localhost/tf/upload.php");
  xhr.send(fd);
	document.getElementById(formId + '-filename').value = xhr.filename;
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

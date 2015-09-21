<?php header('Content-Type: text/javascript; charset=UTF-8'); ?>

$(function(){ // on dom ready
console.info('Start building cytoscape network');
var elesJson = <?php require('createCytoscapeNetwork.php'); echo createJSON('./uploads/' . $_REQUEST['folder'] . '/'.$_REQUEST['dataset'].'/res.txt') .";";?>

console.info(elesJson);
var colors = {

  tf : [ { data : { color : '#CCC' } } ]

};

var elesJsson = {
  nodes: [
    { data: { id: 'a', foo: 3, bar: 5, baz: 7 } },
    { data: { id: 'b', foo: 7, bar: 1, baz: 3 } },
    { data: { id: 'c', foo: 2, bar: 7, baz: 6 } },
    { data: { id: 'd', foo: 9, bar: 5, baz: 2 } },
    { data: { id: 'e', foo: 2, bar: 4, baz: 5 } }
  ], 

  edges: [
    { data: { id: 'ae', weight: 1, source: 'a', target: 'e' } },
    { data: { id: 'ab', weight: 3, source: 'a', target: 'b' } },
    { data: { id: 'be', weight: 4, source: 'b', target: 'e' } },
    { data: { id: 'bc', weight: 5, source: 'b', target: 'c' } },
    { data: { id: 'ce', weight: 6, source: 'c', target: 'e' } },
    { data: { id: 'cd', weight: 2, source: 'c', target: 'd' } },
    { data: { id: 'de', weight: 7, source: 'd', target: 'e' } }
  ]
};


$('#<?php echo $_REQUEST['dataset'];?>' + '-networkPane').cytoscape({
  style: cytoscape.stylesheet()
    .selector('node')
      .css({
        'background-color': 'data(color)',
        'border-width':'1px',
        'border-color': 'black',
        //'mapData(baz, 0, 10, 10, 40)',
        'width': '30px',//'20px',//'mapData(baz, 0, 10, 10, 40)',
        'height': '28px',
        'content': 'data(id)',
        'shape': 'data(shape)'
      })
    .selector('edge')
      .css({
        'line-color': '#888',//'#F2B1BA',
        //'line-color': 'red',
        'target-arrow-color': '#888',//'#F2B1BA',
        'width': 4,
        'target-arrow-shape': 'triangle',
        'opacity': 0.7
      })
    .selector(':selected')
      .css({
        'background-color': 'red',
        'line-color': 'black',
        'target-arrow-color': 'red',
        'source-arrow-color': 'red',
        'opacity': 1
      })
    .selector('.faded')
      .css({
        'opacity': 0.25,
        'text-opacity': 0
      }),
  
  elements: elesJson,
  
  layout: {
    name: 'breadthfirst',
    padding: 10
  },
  
  panningEnabled: true,
  userPanningEnabled:true,
  headless:false,
  styleEnabled:true,
  hideEdgesOnViewport:false,
  hideLabelsOnViewport:false,
  autolock:false,

  ready: function(){
    // ready 1
    console.log("Network for <?php echo $_REQUEST['dataset'];?> has been built");
    this.panzoom();
  }
});


/*
$('#cy2').cytoscape({
  style: cytoscape.stylesheet()
    .selector('node')
      .css({
        'background-color': '#6272A3',
        'shape': 'rectangle',
        'width': 'mapData(foo, 0, 10, 10, 30)',
        'height': 'mapData(bar, 0, 10, 10, 50)',
        'content': 'data(id)'
      })
    .selector('edge')
      .css({
        'width': 'mapData(weight, 0, 10, 3, 9)',
        'line-color': '#B1C1F2',
        'target-arrow-color': '#B1C1F2',
        'target-arrow-shape': 'triangle',
        'opacity': 0.8
      })
    .selector(':selected')
      .css({
        'background-color': 'black',
        'line-color': 'black',
        'target-arrow-color': 'black',
        'source-arrow-color': 'black',
        'opacity': 1
      }),
  
  elements: elesJson,
  
  layout: {
    name: 'breadthfirst',
    directed: true,
    padding: 10
  },
  
  ready: function(){
    // ready 2
  }
});
  */


}); // on dom ready

$(function() {
var <?php echo str_replace('-', '_', $_REQUEST['dataset']); ?>Network = $('#<?php echo $_REQUEST['dataset'];?>' + '-networkPane').cytoscape('get');

$('#<?php echo $_REQUEST['dataset'];?>' + '-networkPane').cytoscape('get').on(
  'select','node', function(evt) {
    var cy = $('.cytoPanel').cytoscape('get');

    var names = "";
    cy.nodes(':selected').forEach(function(ele) {
      names = names + ele.id() + " ";
    });
    $("#selectedNodes").html(names);
  });

$('#<?php echo $_REQUEST['dataset'];?>' + '-networkPane').cytoscape('get').on(
  'unselect','node', function(evt) {
    var cy = $('.cytoPanel').cytoscape('get');

    var names = "";
    cy.nodes(':selected').forEach(function(ele) {
      names = names + ele.id() + " ";
    });
    $("#selectedNodes").html(names);
  });


});
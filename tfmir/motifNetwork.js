$('#motifNetwork').cytoscape({
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
    console.log("Network for has been built");
    //console.log(elesJson);
    //this.panzoom();
    }
  });
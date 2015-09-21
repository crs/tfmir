<?php

$ints = array('Count', 'target.reg', 'regulator.reg','ENTREZ.IDs.coreg.targets','ENTREZ.IDs.cotargeted.genes');
$strings = array('Category','Term','Mir','category',  'regulator','target',	'evidence',	'source',	'regulator.association.disease', 'regulator.association.function', 'target.association.disease', 'target.association.function', 'motifs.ids','motif.type','tf','mirna','targets','Symbol.coreg.targets','Symbol.cotargeted.genes');
$floats = array('Percentage','Pval','Pval.BH','Pval.Bonf','regulator.dsw','target.dsw','pvals.tfmirpair','pval.random','zscore');
$bools = array('is_regulator_in_disease', 'is_target_in_disease');

$types = array('int' => $ints, 'string' => $strings, 'float' => $floats, 'bool' => $bools);

$dic = parseMappings('dict-mappings.txt');

function getKeyType($key) {
  global $types;
  foreach ($types as $type => $typestring) {
    if (in_array($key, $typestring)) { 
      return $type;
    }
  }
  echo 'WARNING: ' .$key . ' did not trigger anything!';
}


function parseMappings($filename) {
  return parse_ini_file($filename);
}

function getValue($key) {
  global $dic;
  if (array_key_exists($key, $dic)) {
    return $dic[$key];
  } else {
    return ucfirst($key);
  }
}

?>
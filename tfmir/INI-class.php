<?php

/*  Copyright 2011 Jason Marshall

This file is part of INI-class.

INI-Class is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

INI-Class is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with INI-Class.  If not, see <http://www.gnu.org/licenses/>.

*/
error_reporting(E_ERROR | E_WARNING | E_PARSE);

define('INI_AUTOSAVE',1); // 0 is not to save at the close of the class, 1 saves
//I really don't know of anywhere else to have an option.  but this is where I would put it.

class parseINI
{
	
	public $ini;		// The INI string, as a full string.  Comments and whitespace stripped.
	public $ini_file;	// The filename.
	public $ini_hand; // The INI file handle.  
	public $ini_array;// The INI broken down into an array
	
	private $ini_opt; // Options array
	
	function __construct($file)
	{
		$this->iniopt();
		$this->ini_file=$file;
		if(!is_file($file)) 
		{
			$brad=fopen($file,'w');
			fwrite($brad,' ');
			fclose($brad);
		} else {
			$inistring=file_get_contents($file);
			$this->ini=trim(preg_replace('/(\s*[#;].*$)|(^\s*\n)/m', '', $inistring));
			$x=$this->parseString();
		}
	}
	
	function parseString($string = null)
	{
		$this->ini_array = array();
		if($string == null) { $string=$this->ini; }
		$jones=preg_split('/\n/is',$string);
		//print_r($jones);
		foreach ($jones as $txt)
		{
			$c=preg_match("/(\\[.*?\\])/is", $txt,$matches);
			if($c)
			{
				$foo=preg_replace($pattern= "/(\\[)/is", $replacement = '', $txt);
				$section=preg_replace($pattern= "/(\\])/is", $replacement = '', $foo); 
				if(!isset($this->ini_array[$section])) { $this->ini_array[$section] = array(); }  
			} else {
				$janet=explode($delimiter='=', $txt, $limit = 2); // My Easter Egg.  I love Rocky Horror
				$colombia=trim($janet[0]);
				$riff=trim($janet[1]);
				$this->ini_array[$section][$colombia]=$riff;
			}  	
		}
		return $this->ini_array;  
	}
	
	
	function parseArray($array = null)
	{
		if ($array == null) { $array=$this->ini_array; }
		$str="#INI Files Created by INI-Class.\n"; //Just something to set the string file up.
		foreach($array as $section => $nxt_arry)
		{  $str=$str."[".$section."]\n";
			foreach($nxt_arry as $name => $value)
			{  $str=$str.$name."=".$value."\n"; }
		}
		$this->ini=$str;
		return($str);
	}
	
	function saveINI($str = null)
	{
		if($str == null) { $this->parseArray(); $str=$this->ini; }
		file_put_contents($this->ini_file,$str);
	}
	
	function iniopt($opt=null,$value=null)
	{
		if($opt==null)
		{
			$this->ini_opt[1]=1;
		} else {
			$this->ini_opt[$opt]=$value;
		}
	}
	
	
	function __destruct()
	{
		if ($this->ini_opt[1]==1) { $this->saveINI($this->ini); }
	}
	
}
?>
<?php
require __DIR__ . "/functions.php";
$type = $argv[1];
$cutoffs = require "cutoffs.php";

$cutoff = $cutoffs[$type];
$files = glob("output.bed.*");
$results = [];
foreach ($files as $file) {
    $feature = features($file);
    foreach ($feature as $f) {
        $results["{$f['chr']}\t{$f['start']}\t{$f['end']}"]++;
    }
}
arsort($results);
$results = array_filter($results, function ($v) use ($cutoff) {
    return $v >= $cutoff;
});
$fp = STDOUT;
foreach ($results as $r => $v) {
    // use fputcsv to escape special characters
    fputcsv($fp, explode("\t", $r), "\t");
}

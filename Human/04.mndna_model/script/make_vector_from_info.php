<?php
require "functions.php";
$info_file = $argv[1];
$mode_name = $argv[2];
$whitelist_file = $argv[3];
$output_name = $argv[4];

$lines = file($info_file, FILE_IGNORE_NEW_LINES);
$total = $lines[0];
$ids = $pos_ids = $neg_ids = [];
for ($i = 1; $i <= $total; $i++) {
    $line = $lines[$i];
    $line = explode(" ", $line);
    $id = $line[0];
    $label = $line[1];
    if ($label == 1) {
        $pos_ids[] = $id;
    } else {
        $neg_ids[] = $id;
    }
    $ids[] = $id;
}
$whitelist = file($whitelist_file, FILE_IGNORE_NEW_LINES);
natsort($whitelist);
$dir = dirname($output_name);
$window_selection = "$dir/learn.window.selection";
file_put_contents($window_selection, count($whitelist) . " 10000" . "\n" . join("\n", $whitelist) . PHP_EOL);
make_vector($ids, $mode_name, $output_name, [], $whitelist, $pos_ids, $neg_ids);

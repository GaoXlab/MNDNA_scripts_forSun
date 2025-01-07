<?php

require __DIR__ . "/functions.php";

$train_file_name = $argv[1];
$test_file_name = $argv[2];
$mode = $argv[3];
$whitelist_file = $argv[4];
$test_label = [];

$train_ids = load_id_from_file($train_file_name, $train_label);
$test_ids = load_id_from_file($test_file_name, $test_label);
$train_ids = array_diff($train_ids, $test_ids);

$whitelist = file($whitelist_file, FILE_IGNORE_NEW_LINES);
// always sort whitelist by natural order
natsort($whitelist);
file_put_contents("learn.window.selection", count($whitelist) . " " . 10000 . PHP_EOL);
file_put_contents("learn.window.selection", join("\n", $whitelist), FILE_APPEND);
$samples = [];
foreach ($train_label as $id => $label) {
    $samples[$id] = ['ca' => $label];
}
foreach ($test_label as $id => $label) {
    $samples[$id] = ['ca' => $label];
}
make_vector($train_ids, $mode, "learn.vector.train", $samples, $whitelist);
make_vector($test_ids, $mode, "learn.vector.test", $samples, $whitelist);
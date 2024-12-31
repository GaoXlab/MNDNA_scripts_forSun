<?php

require __DIR__ . "/functions.php";

$pos_file_name = $argv[1];
$neg_file_name = $argv[2];
$test_file_name = $argv[3];
$mode = $argv[4];
$whitelist_file = $argv[5];
$pos_ids = load_id_from_file($pos_file_name);
$neg_ids = load_id_from_file($neg_file_name);
$test_label = [];
$test_ids = load_id_from_file($test_file_name, $test_label);
$whitelist = file($whitelist_file, FILE_IGNORE_NEW_LINES);
// always sort whitelist by natural order
natsort($whitelist);
file_put_contents("learn.window.selection", count($whitelist) . " " . 10000 . PHP_EOL);
file_put_contents("learn.window.selection", join("\n", $whitelist), FILE_APPEND);
$train_ids = array_merge($pos_ids, $neg_ids);
$train_ids = array_diff($train_ids, $test_ids);
$samples = [];
$test_pos_ids = array_keys(array_filter($test_label, function ($v) {
    return $v == 1;
}));
$test_neg_ids = array_keys(array_filter($test_label, function ($v) {
    return $v == 0;
}));

make_vector($train_ids, $mode, "learn.vector.train", $samples, $whitelist, $pos_ids, $neg_ids);
make_vector($test_ids, $mode, "learn.vector.test", $samples, $whitelist, $test_pos_ids, $test_neg_ids);
<?php
require __DIR__ . "/functions.php";
$files = glob("output.bed.*");
$results = [];
$type = $argv[1];
$ids_file = $argv[2];

$cutoffs = [
    'crc' => ['width' => 30, 'median' => 30],
    'panca' => ['width' => 10, 'median' => 10],
];
foreach ($files as $file) {
    $feature = features($file);
    foreach ($feature as $f) {
        $results["{$f['chr']}\t{$f['start']}\t{$f['end']}"]++;
    }
}
$cutoff = $cutoffs[$type];
arsort($results);
$feature_detail = get_feature_detail($type, array_keys($results), $ids_file);
$results = array_filter($results, function ($v, $k) use ($cutoff, $feature_detail) {
    return $v >= 10 && $feature_detail[$k]['width'] > $cutoff['width'] && $feature_detail[$k]['median'] >= $cutoff['median'];
}, ARRAY_FILTER_USE_BOTH);
$fp = STDOUT;
foreach ($results as $r => $v) {
    // use fputcsv to escape special characters
    fputcsv($fp, explode("\t", $r), "\t");
}

function get_feature_detail($type, $features, $ids_file)
{
    $tab_file_dir = dirname($ids_file);
    foreach($features as $feature) {
        list($chr, $start, $end) = explode("\t", $feature);
        $feature_detail[$feature] = [
            'width' => ($end - $start)/1000,
        ];
    }
    $ids = load_id_from_file($ids_file);
    $index = file("$tab_file_dir/manu_{$type}_240103/sorted.tab.index", FILE_IGNORE_NEW_LINES);
    foreach ($ids as $id) {
        $data[$id] = array_combine($index, raw($id, "manu_{$type}_240103"));
    }
    foreach ($feature_detail as $feature => $detail) {
        $feature_detail[$feature]['median'] = median(array_column($data, $feature));
    }
    return $feature_detail;
}

function median(array $array_column)
{
    sort($array_column);
    $count = count($array_column);
    $middle = floor(($count - 1) / 2);
    if ($count % 2) {
        return $array_column[$middle];
    } else {
        return ($array_column[$middle] + $array_column[$middle + 1]) / 2;
    }
}
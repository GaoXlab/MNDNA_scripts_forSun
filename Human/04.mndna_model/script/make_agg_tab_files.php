<?php
$agg = $argv[1];
$black_list_file = __DIR__ . "/blacklist.bed";

$input_file_name = "train.tab";
$output_filename = "{$input_file_name}.{$agg}";

$fp = fopen($input_file_name, 'r');

$head = fgetcsv($fp, null, "\t", "'");
$seq_id = array_map(function ($v) {
    return explode('.', $v)[0];
}, array_slice($head, 3));
fclose($fp);
$fp = fopen($input_file_name, 'r');
$fp_output = fopen($output_filename, 'w+');

# copy input header to output file
$head = fgetcsv($fp, null, "\t");
fputcsv($fp_output, $head, "\t");

$col = count($head);

$pos = [];
$sum = [];
$sum[0] = array_fill(0, count($head), 0);
$i = 1;
while ($line = fgetcsv($fp, null, "\t")) {
    $pos[$i] = array_slice($line, 0, 3);
    $pos[$i][2] = $pos[$i][1] + $agg * 10000;
    for ($j = 3; $j < $col; $j++) {
        $sum[$i][$j] = $sum[$i - 1][$j] + $line[$j];
    }
    if ($i >= $agg) {
        $output_index = $i - $agg + 1;
        $output_data = [];
        for ($j = 3; $j < $col; $j++) {
            $output_data[] = $sum[$i][$j] - $sum[$output_index - 1][$j];
        }
        unset($sum[$output_index - 1]);
        fputcsv($fp_output, [...$pos[$output_index], ...$output_data], "\t");
    }
    $i++;
    if ($i % 50000 == 1) {
        echo "$output_filename $i \n";
    }
}
exec("head -n 1 $output_filename > $output_filename.clean && bedtools subtract -a $output_filename -b $black_list_file -A >> $output_filename.clean && mv $output_filename.clean $output_filename");
echo "$output_filename DONE\n";

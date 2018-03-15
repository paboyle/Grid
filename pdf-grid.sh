echo 'Kill all Jekyll instances'
kill -9 $(ps aux | grep '[j]ekyll' | awk '{print $2}')
clear

echo "Building PDF-friendly HTML site for Grid ...";
bundle exec jekyll serve --detach --config _config.yml,pdfconfigs/config_grid_pdf.yml;
echo "done";

sed -i "/^\s*$/d" _site/pdfconfigs/prince-list.txt 
sed 's/http:\/\/localhost:4010\/grid-pdf/_site/' _site/pdfconfigs/prince-list.txt > site_list

echo "Preprocess generated html files..."
cat site_list | while read filename;
do
    complete_filename=$filename
    if ! [[ $filename = *".html" ]]; then
	complete_filename="${filename}index.html"
    fi
    echo $complete_filename
    phantomjs assets/js/render-math-and-dump.js $complete_filename | sed -n '/DOCTYPE/,$p' > temp 
    mv temp $complete_filename
done
rm site_list
echo "done"

echo "Building the PDF ...";
prince --javascript --no-warn-css --input-list=_site/pdfconfigs/prince-list.txt -o pdf/grid.pdf;

echo "Done. Look in the pdf directory to see if it printed successfully."

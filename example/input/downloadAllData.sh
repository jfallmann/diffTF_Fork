# This script can be used to download all required data for the example analysis.
# Simply execute with "sh downloadAllData.sh"

filename="exampleData.tar.gz"
wget -O $filename https://www.embl.de/download/zaugg/diffTF/$filename && tar xvzf $filename --overwrite && rm $filename

filename="mm10.fa.tar.gz"
wget -O $filename https://www.embl.de/download/zaugg/diffTF/referenceGenome/$filename && mkdir referenceGenome && tar xvzf $filename -C referenceGenome --overwrite && rm $filename

filename="TFBS_mm10_PWMScan_HOCOMOCOv10.tar.gz"
wget -O $filename https://www.embl.de/download/zaugg/diffTF/TFBS/$filename && tar xvzf $filename --overwrite && rm $filename

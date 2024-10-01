vcfStripDttm(){
    local path=$1

    zgrep -v "^##fileDate" ${path} | gzip -n -f > tmp.gz && mv tmp.gz ${path}

}


gzipStripDttm(){
    local path=$1

    mv ${path} tmp.gz 
    gunzip -f tmp.gz && gzip -n -f tmp 
    mv tmp.gz ${path}

}

extractDemuxTsv(){
    local infile=$1
    local library=$2

    awk {'printf ("%s\t%s\t%s\t%s\t%s\n", $2, $3, $4, $5, $6)'} ${infile} > ${library}.clust1.samples.reduced.tsv

}

demuxTsvFromFmx(){
    local library=$1

    gunzip -f ${library}.clust1.samples.gz
    extractDemuxTsv ${library}.clust1.samples ${library}
    gzip -f -n ${library}.clust1.samples
}
for L in {A..Z}; do
 for list in `ls -1 ~/pgba_work/ouleye/PPR_MGEproject/migale_data_1/$L*/*`; do
   sp=`echo $list| cut -d '/' -f 8`
   gcf=`echo $list| cut -d '/' -f '9'`
   echo $list >> all_paths.txt
 
   gzip -d -c /db/gb_bacteria/current/flat/$sp/latest_assembly_versions/$gcf/${gcf}_wgsmaster.gbff.gz| awk 'BEGIN{src="";host="";FS="\""}{if ($0 ~ "isolation_source"){src=$2};if ($0 ~ "/host="){host=$2}}END{print src"\t"host}'
 
  done
done

#! /bin/gawk -f

BEGIN {FS="\t"}

# Sample categorization
{
 if (tolower($0) ~ "gut|intestin|fece|faece|fecal|faecal|feace|stool|caecum|caeca|ceacal|ceca|cecum|duodenum|jejunum|rumen|colon|digesti|ileum|ileal|rectum|rectal|cloaca|chyme|manure")
   scat="gut"
 else
 if ($1 =="" || tolower($1) ~ "not provided|not available|biological product|rich broth|glycerol stored|clinical isolate|anaerobe media|^swab$|^tissue$"|| tolower($1) ~ "^(chicken|hen|broiler|fowl|poultry|rooster|human|homo sapiens|man|woman|infant|duck|goose|mouse|rat|pig|swine|sow|pork|cow|cattle|calf|rabbit|bunny|hamster|guinea pig|cavy|mare|horse|poney|donkey|mule|equine|cat|dog)$")
  scat="unknown"
 else
  scat="non gut"

# Host categorization
 if (tolower($0) ~ "honey bee|pigeon|pheasant")
   hcat="other"
 else  
 if (tolower($0) ~ "duck|(^|\\W|$)Anas(^|\\W|$)|turkey|meleagris|goose|(^|\\W|$)anser(^|\\W|$)|mouse|mus musculus|sheep|(^|\\W|$)lamb(^|\\W|$)|(^|\\W|$)ovis(^|\\W|$)|goat|(^|\\W|$)capra(^|\\W|$)|(^|\\W|$)rat(^|\\W|$)|rattus|rabbit|oryctolagus|bunny|hamster|mesocricetus|guinea pig|cavy|(^|\\W|$)cavia(^|\\W|$)|quail|(^|\\W|$)mare(^|\\W|$)|horse|poney|donkey|(^|\\W|$)mule(^|\\W|$)|equine|equus|(^|\\W|$)cat(^|\\W|$)|felis|(^|\\W|$)dog(^|\\W|$)|canis lupus familiaris|alpaca|vicugna pacos")
   hcat="other domestic"
 else
 if (tolower($0) ~ "(^|\\W|$)pig|swine|sus scrofa|suidae|(^|\\W|$)hog(^|\\W|$)|sow|pork")
   hcat="pig"
 else  
 if (tolower($0) ~ "cow|cattle|bos taurus|bovin|heifer|calf|beef|steer|charolais|holstein|limousine|angus|hereford")
   hcat="cattle"
 else  
 if (tolower($0) ~ "gallus gallus|chick|(^|\\W|$)hen(^|\\W|$)|broiler|fowl|cockerel|poultry|rooster")
   hcat="chicken"
 else
  if (tolower($0) ~ "human|homo sapiens|hoimo sapiens|(^|\\W|$)man(^|\\W|$)|woman|infant|adult|baby|people")
   hcat="human"
 else
 if ($2 =="" || tolower($2) ~ "not provided|^animal$")
   hcat="unknown"
 else
   hcat="other"

   
   
 print $0"\t"scat"\t"hcat  
}
   


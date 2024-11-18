#! /usr/bin/gawk -f

#gawk -f EssaipourGhislain.awk -v path="/myfolder/" -v nosimu="8" -v nochr="1,2" -v phasing="true" oldsettings.toml


BEGIN{
}
{
if($1=="dir"){
    gsub(/dir_path/, path, $3)
    {print $1," ",$2," ",$3}
} else if($1=="file" || $1=="Nea" || $1=="CEU" || $1=="YRI") {
    gsub(/dir_path/, path, $3)
    gsub(/simid/, nosimu, $3)
    {print $1," ",$2," ",$3}
} else if($1=="chr") {
    gsub(/dir_path/, path, $3)
    gsub(/chrlist/, nochr, $3)
    {print $1," ",$2," ",$3}
} else if($1=="phased") {
    gsub(/hap_phasing/, phasing, $3)
    {print $1," ",$2," ",$3}
} else {print}

}END{
}
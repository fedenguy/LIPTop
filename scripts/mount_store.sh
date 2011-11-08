base_storedir="./store"
storeremote="pclip11.cern.ch:/localdata/`whoami`/"
storedir="${base_storedir}"
echo $storeremote
if ( ! mount | grep ${storeremote} > /dev/null )
then 
    [[ ! -d ${storedir} ]] && mkdir ${storedir}
    sshfs ${storeremote} ${storedir} -o nonempty
fi



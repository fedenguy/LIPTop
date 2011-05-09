base_storedir="${CMSSW_BASE}/datastore"
storedir="${base_storedir}"
storeremote="ttbar2l@pclip11.cern.ch:/localdata/"

if ( ! mount | grep ${storeremote} > /dev/null )
then 
    [[ ! -d ${storedir} ]] && mkdir ${storedir}
    sshfs ${storeremote} ${storedir} -o nonempty
fi



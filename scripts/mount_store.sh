base_storedir="./store"
#storeremote="pclip11.cern.ch:/localdata/`whoami`/"
storeremote="pclip11.cern.ch:/localdata/"
WHOAMI=`whoami`
if [ "$WHOAMI"="psilva" ];
    then
    storeremote="maccms316.cern.ch:/Users/psilva/store/"
#    storeremote="cmsphys03.cern.ch:/data/psilva/"
fi
echo "Mounting ${storeremote} for ${WHOAMI} @ ${base_storedir}"

storedir="${base_storedir}"
echo $storeremote
if ( ! mount | grep ${storeremote} > /dev/null )
then 
    [[ ! -d ${storedir} ]] && mkdir ${storedir}
    sshfs ${storeremote} ${storedir} -o nonempty
else
    echo "Failed to mount"
fi



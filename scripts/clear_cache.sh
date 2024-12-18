

#!/bin/bash

# Define the password variable
SUDO_PASSWORD='David_24'

# Function to clear cache
clear_cache() {
    echo $SUDO_PASSWORD | sudo -S sync
    echo $SUDO_PASSWORD | sudo -S sysctl -w vm.drop_caches=3
    echo "Cache cleared at $(date)"
}

# Run the cache clear function every 30 minutes
while true; do
    clear_cache
    sleep 600  # Sleep for 30 minutes (1800 seconds)
done

#sudo sh -c "sync; echo 3 > /proc/sys/vm/drop_caches"

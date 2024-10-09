#!/bin/bash

set -x

# We start by adding extra apt packages, since pip modules may required library
if [ "$EXTRA_APT_PACKAGES" ]; then
    echo "EXTRA_APT_PACKAGES environment variable found.  Installing."
    apt update -y
    apt install -y $EXTRA_APT_PACKAGES
fi

if [ "$USE_MAMBA" == "true" ]; then
    echo "USE_MAMBA enabled. Using mamba for all conda operations"
    CONDA_BIN="/opt/conda/bin/mamba"
else
    CONDA_BIN="/opt/conda/bin/conda"
fi

if [ -e "/opt/app/environment.yml" ]; then
    echo "environment.yml found. Installing packages"
    $CONDA_BIN env update -f /opt/app/environment.yml
else
    echo "no environment.yml"
fi

if [ "$EXTRA_CONDA_PACKAGES" ]; then
    echo "EXTRA_CONDA_PACKAGES environment variable found.  Installing."
    $CONDA_BIN install -y $EXTRA_CONDA_PACKAGES
fi

if [ "$EXTRA_PIP_PACKAGES" ]; then
    echo "EXTRA_PIP_PACKAGES environment variable found.  Installing".
    /opt/conda/bin/pip install $EXTRA_PIP_PACKAGES
fi



echo "CLUSTER_NAME=$CLUSTER_NAME"
echo "NAMESPACE=$NAMESPACE"
SCHEDULER_ADDRESS="${CLUSTER_NAME}-dask-scheduler.${NAMESPACE}.svc.cluster.local"
echo "SCHEDULER_ADDRESS=$SCHEDULER_ADDRESS"
echo "SCHEDULER_PORT=$SCHEDULER_PORT"
echo "SERVER_PORT=$SERVER_PORT"
echo "STORAGE_PATH=$STORAGE_PATH"

conda list

mkdir -p $STORAGE_PATH/working-directory
/opt/conda/bin/python /usr/bin/server.py --scheduler-address $SCHEDULER_ADDRESS --scheduler-port $SCHEDULER_PORT --cluster-name $CLUSTER_NAME --port $SERVER_PORT --storage-path $STORAGE_PATH/working-directory

# Run extra commands
exec "$@"
# Running Evaluator on NRP with Helm

This is a proof-of-concept execution of Evaluator on NRP.

## Set up persistent storage

```
kubectl apply -f pvc.yaml
```

Check it's bound:

```
kubectl get pvc -o wide
```


## Deploy a Dask cluster using Helm

Customize `helm-values.yaml` as you wish and run the below. You will need to be an admin of your namespace:

```
MY_CLUSTER_NAME="evaluator-lw"
helm install -f helm-values.yaml --repo https://helm.dask.org  $MY_CLUSTER_NAME dask
```

If you get an error like `cannot re-use a name that is still in use`:

```
helm install -f helm-values.yaml --repo https://helm.dask.org --replace $MY_CLUSTER_NAME dask

```

You will get output similar to the below. Execute the port-forwarding commands.

```
NAME: evaluator-lw
LAST DEPLOYED: Fri Oct  4 16:30:04 2024
NAMESPACE: openforcefield
STATUS: deployed
REVISION: 1
TEST SUITE: None
NOTES:
Thank you for installing DASK, released at name: evaluator-lw.

To learn more about the release, try:

  $ helm status evaluator-lw  # information about running pods and this message
  $ helm get all evaluator-lw     # get full Kubernetes specification

This release includes a Dask scheduler, 3 Dask workers, and  Jupyter servers.

The Jupyter notebook server and Dask scheduler expose external services to
which you can connect to manage notebooks, or connect directly to the Dask
cluster. You can get these addresses by running the following:

  export DASK_SCHEDULER="127.0.0.1"
  export DASK_SCHEDULER_UI_IP="127.0.0.1"
  export DASK_SCHEDULER_PORT=8080
  export DASK_SCHEDULER_UI_PORT=8081
  kubectl port-forward --namespace openforcefield svc/evaluator-lw-dask-scheduler $DASK_SCHEDULER_PORT:8786 &
  kubectl port-forward --namespace openforcefield svc/evaluator-lw-dask-scheduler $DASK_SCHEDULER_UI_PORT:80 &

  export JUPYTER_NOTEBOOK_IP="127.0.0.1"
  export JUPYTER_NOTEBOOK_PORT=8082
  kubectl port-forward --namespace openforcefield svc/evaluator-lw-dask-jupyter $JUPYTER_NOTEBOOK_PORT:80 &

  echo tcp://$DASK_SCHEDULER:$DASK_SCHEDULER_PORT               -- Dask Client connection
  echo http://$DASK_SCHEDULER_UI_IP:$DASK_SCHEDULER_UI_PORT     -- Dask dashboard
  echo http://$JUPYTER_NOTEBOOK_IP:$JUPYTER_NOTEBOOK_PORT       -- Jupyter notebook

NOTE: It may take a few minutes for the LoadBalancer IP to be available. Until then, the commands above will not work for the LoadBalancer service type.
You can watch the status by running 'kubectl get svc --namespace openforcefield -w evaluator-lw-dask-scheduler'
```

## Set up the Evaluator Server

Customize the pod as you choose and set up the Evaluator server as a separate pod.

```
kubectl apply -f server-pod.yaml
```

## Running

Forward whichever port you choose for the Evaluator server (default 8998).

```
kubectl port-forward --namespace openforcefield pod/evaluator-server-pod 8998:8998 &
```

Run the script.
```
python run-client.py
```

### Debugging

To download working files:

```
kubectl cp evaluator-server-pod:/evaluator-storage/working-directory .
```


### Notes for scaling up (or down)

Since May 2024, Dask-Kubernetes has removed the HelmCluster that allowed scaling within Python.
To manually scale with Kubernetes:

```
kubectl scale --replicas=1 deployment/${MY_CLUSTER_NAME}-dask-worker
```


## Deleting

Stop the server:

```
kubectl delete -f server-pod.yaml
```

Remove storage:

```
kubectl delete -f pvc.yaml
```

Stop Dask:

```
helm delete $MY_CLUSTER_NAME
```


If you want your ports back, kill the processes where kubectl is forwarding them.

```
lsof -ti:8998 | xargs kill -9
```

(I have the below in a `.bash_aliases` file that is sourced by `.bashrc`; you might find it useful.)

```
function killport {
    port=$1
    lsof -ti:${port} | xargs kill -9
}
```



## How it works

The local EvaluatorClient passes data to the remote EvaluatorServer, communicating through the forwarded port.
The EvaluatorServer communicates with the Dask scheduler launched via Helm to handle compute.
The persistent volume claim is necessary as currently data is stored via files passed via the filesystem to future protocols. Both the EvaluatorServer and dask workers must be connected to the same file system.
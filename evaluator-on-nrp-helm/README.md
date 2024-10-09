# Running Evaluator on NRP with Helm

This is a proof-of-concept execution of Evaluator on NRP.

## Set up persistent storage

We will need a persistent file system. Set one up below:

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
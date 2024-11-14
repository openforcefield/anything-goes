# Run Evaluator on KubeCluster

The files in this directory allow a proof-of-concept Evaluator fit to run on NRP.

Everything happens in `run-job.py`. Run as below:

```
python run-job.py > log.log 2>&1
```

It does the following:

* Creates a PVC for persistent file storage across workers
* Create an calculation backend using a [Dask KubeCluster](https://kubernetes.dask.org/en/latest/operator_kubecluster.html) as a backend. (It writes out `cluster-spec.yaml` if you want to have a look at what's being applied.)
* Copy the `server-existing.py` script across
* Create a Kubernetes deployment running an EvaluatorServer. The EvaluatorServer *connects directly to the Dask KubeCluster* through the local TCP address. Unfortunately, we don't have sufficient permissions to connect through the Dask KubeCluster API.
* Forward the Evaluator server port locally to the port specified in `targets/phys-prop/options.json`
* Run `ForceBalance.py optimize.in`
* On an exception or interruption, cleans up the EvaluatorServer Deployment and the PVC. The DaskCluster should exit itself on an error or interruption.

**Note:**

The options for connecting to the EvaluatorServer are in `targets/phys-prop/options.json` and are used through the ForceBalance Evaluator interface. For this experimental fit, if you want to change the port from 8998 to something else, you will need to edit that file too.

## Keeping an eye on progress

If you want to have a look at progress in the working directory and track how large files are getting, create a Pod to interface with the storage. The `copy-data.sh` script copies the working directory locally so you can check on progress.

```
kubectl apply -f pod-interface.yaml
./copy-data.sh
```


## Cleaning up manually

The try/finally only... sometimes works.
In particular, it may fail if you've created additional pods to interact with the PVC.
If you need to clean up the kubernetes resources manually, look for:

```
kubectl get daskclusters
kubectl get deployments
kubectl get pvc
```

And delete:

```
kubectl delete daskcluster evaluator-lw
kubectl delete deployment evaluator-server-lw
kubectl delete pvc evaluator-storage-lw
```

Note that it's important for the PVC deletion to happen last.

Also, you may need to manually stop port-forwarding your local port.



## Things I tried that didn't work

### RBAC

A *lot* of this convoluted code comes from the need for the EvaluatorServer to run remotely on NRP. Dask-Kubernetes and associated libraries are not well-set-up for users who do not have cluster-wide permissions. You will get an error like:

```
User [xxx] cannot get resource "clusterrolebindings" in API group "rbac.authorization.k8s.io" at the cluster scope
```

if you try to follow the Dask RBAC guide.

### kr8s

kr8s is a much nicer library to use than straight kubernetes, but unfortunately also ran into permissions issues when authenticating:

```
User [xxx] cannot create resource "tokenreviews" in API group "authentication.k8s.io" at the cluster scope
```

This is a bit mysterious since it looks like modern Dask kubernetes does use kr8s under the hood. I might want to look into this in the future.

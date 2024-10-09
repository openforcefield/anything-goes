import argparse
import time
import logging

from openff.evaluator.backends.dask import BaseDaskBackend, _Multiprocessor, QueueWorkerResources
from openff.evaluator.server import EvaluatorServer
import distributed

logger = logging.getLogger(__name__)

parser = argparse.ArgumentParser(description="Run Evaluator server")
parser.add_argument(
    "--scheduler-address", "-sa",
    type=str,
    default="evaluator-lw-dask-scheduler.openforcefield.svc.cluster.local",
)
parser.add_argument(
    "--scheduler-port", "-sp",
    type=int,
    default=8786,
)
parser.add_argument(
    "--cluster-name", "-n",
    type=str,
    default="evaluator-lw",
)
parser.add_argument(
    "--port", "-p",
    type=int,
    default=8998,
)
parser.add_argument(
    "--storage-path",
    type=str,
    default="/evaluator-storage/working-directory"
)

class DaskHelmInterface(BaseDaskBackend):
    def __init__(
        self,
        
        address: str = "evaluator-lw-dask-scheduler.openforcefield.svc.cluster.local",
        port: int = 8786,  # $DASK_SCHEDULER_PORT
        name: str = "evaluator-lw",
    ):
        number_of_workers = 2
        resources_per_worker = QueueWorkerResources()
        super().__init__(number_of_workers, resources_per_worker)

        assert isinstance(resources_per_worker, QueueWorkerResources)

        if resources_per_worker.number_of_gpus > 0:
            if resources_per_worker.number_of_gpus > 1:
                raise ValueError("Only one GPU per worker is currently supported.")
            
        self._cluster = f"tcp://{address}:{port}"
        self._name = name

    @staticmethod
    def _wrapped_function(function, *args, **kwargs):
        available_resources = kwargs["available_resources"]
        gpu_assignments = kwargs.pop("gpu_assignments")

        if available_resources.number_of_gpus > 0:
            worker_id = distributed.get_worker().id

            available_resources._gpu_device_indices = (
                "0" if worker_id not in gpu_assignments else gpu_assignments[worker_id]
            )

            logger.info(
                f"Launching a job with access to GPUs "
                f"{available_resources._gpu_device_indices}"
            )

        return_value = _Multiprocessor.run(function, *args, **kwargs)
        return return_value
    
    def stop(self):
        logger.warning(
            f"You will need to stop the cluster manually with ``helm delete {self._name}``"
        )

    def submit_task(self, function, *args, **kwargs):
        key = kwargs.pop("key", None)

        return self._client.submit(
            self._wrapped_function,
            function,
            *args,
            **kwargs,
            key=key,
            available_resources=self._resources_per_worker,
            gpu_assignments={},
        )
    

if __name__ == "__main__":
    args = parser.parse_args()
    print(args)
    calculation_backend = DaskHelmInterface(
        address=args.scheduler_address,
        port=args.scheduler_port,
        name=args.cluster_name,
    )

    print(calculation_backend)
    logger.info(f"Calculating with backend {calculation_backend}")
    logger.info(f"Starting server on port {args.port}")
    with calculation_backend:
        evaluator_server = EvaluatorServer(
            calculation_backend,
            working_directory=args.storage_path,
            port=args.port,
            delete_working_files=False,
        )
        logger.info(f"Server started on port {args.port}")
        evaluator_server.start(asynchronous=False)
        logger.info(f"Server stopped on port {args.port}")

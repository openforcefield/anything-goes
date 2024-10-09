
import argparse
import logging
import json

import distributed

from openff.evaluator.datasets import PhysicalPropertyDataSet
from openff.evaluator.client import RequestOptions
from openff.evaluator.client import EvaluatorClient, ConnectionOptions
from openff.evaluator.forcefield import SmirnoffForceFieldSource
from openff.evaluator.properties import Density

from openff.evaluator.server import EvaluatorServer
from openff.evaluator.backends.dask import BaseDaskBackend, _Multiprocessor
from openff.evaluator.server import EvaluatorServer
from openff.evaluator.backends.backends import ComputeResources, QueueWorkerResources

logger = logging.getLogger(__name__)


parser = argparse.ArgumentParser(description="Run Evaluator server")
parser.add_argument(
    "--scheduler-address", "-sa",
    type=str,
    default="127.0.0.1",
)
parser.add_argument(
    "--scheduler-port", "-sp",
    type=int,
    default=8080,
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

class DaskHelmInterface(BaseDaskBackend):
    def __init__(
        self,
        address: str = "127.0.0.1",
        port: int = 8080,  # $DASK_SCHEDULER_PORT
        name: str = "evaluator-lw",
        resources_per_worker=QueueWorkerResources(
            number_of_gpus=1,
            preferred_gpu_toolkit=ComputeResources.GPUToolkit.CUDA,
        ),
    ):
        super().__init__(2, resources_per_worker)

        assert isinstance(resources_per_worker, QueueWorkerResources)

        if resources_per_worker.number_of_gpus > 0:
            if resources_per_worker.number_of_gpus > 1:
                raise ValueError("Only one GPU per worker is currently supported.")

        self._cluster = f"tcp://{address}:{port}"
        self._name = name
        self._started = False

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

def main(
    port: int = 8000,
    scheduler_address: str = "127.0.0.1",
    scheduler_port: int = 8080,
    cluster_name: str = "evaluator-lw",

    data_set_path: str = "small-test-set.json",
    force_field_path: str = "openff-2.2.0.offxml",

    working_directory: str = "/evaluator-storage/working-directory",
    enable_data_caching: bool = True,
):
    import os

    data_set = PhysicalPropertyDataSet.from_json(data_set_path)
    force_field_source = SmirnoffForceFieldSource.from_path(force_field_path)

    density_schema = Density.default_simulation_schema(n_molecules=256)
    for schema in density_schema.workflow_schema.protocol_schemas:
        # equilibration_simulation is a Protocol                                                                                                                                
        if "simulation" in schema.id:
            schema.inputs[".steps_per_iteration"] = 10
            schema.inputs[".output_frequency"] = 10
        # conditional_group is a group of protocols                                                                                                                             
        if "conditional" in schema.id:
            for protocol_name, protocol in schema.protocol_schemas.items():
                if "simulation" in protocol_name:
                    protocol.inputs[".steps_per_iteration"] = 10
                    protocol.inputs[".output_frequency"] = 10
    
    # Create an options object which defines how the data set should be estimated.                                                                                              
    estimation_options = RequestOptions()
    # Specify that we only wish to use molecular simulation to estimate the data set.                                                                                           
    estimation_options.calculation_layers = ["SimulationLayer"]

    # Add our custom schemas, specifying that the should be used by the 'SimulationLayer'                                                                                       
    estimation_options.add_schema("SimulationLayer", "Density", density_schema)



    # backend = DaskHelmInterface(
    #     port=scheduler_port,
    #     address=scheduler_address,
    #     name=cluster_name,
    # )

    # with backend:
    #     server = EvaluatorServer(
    #         calculation_backend=backend,
    #         working_directory=working_directory,
    #         port=port,
    #         enable_data_caching=enable_data_caching,
    #     )
    #     with server:
    connection_options = ConnectionOptions(
        server_address="localhost",
        server_port=port,
    )
    evaluator_client = EvaluatorClient(
        connection_options=connection_options
    )

    request, exception = evaluator_client.request_estimate(
        property_set=data_set,
        force_field_source=force_field_source,
        options=estimation_options,
    )
    assert exception is None, exception

    results, exception = request.results(synchronous=True, polling_interval=30)
    assert exception is None, exception
    print(f"# queued: {len(results.queued_properties)}")

    print(f"# estimated: {len(results.estimated_properties)}")

    print(f"# unsuccessful: {len(results.unsuccessful_properties)}")
    print(f"# exceptions: {len(results.exceptions)}")
    results.estimated_properties.json("estimated_data_set.json", format=True)
    print(results.exceptions)
    with open("exceptions.json", "w") as f:
        json.dump(results.exceptions, f, indent=2)



if __name__ == "__main__":
    args = parser.parse_args()
    main(
        port=args.port,
        scheduler_address=args.scheduler_address,
        scheduler_port=args.scheduler_port,
        cluster_name=args.cluster_name,
    )

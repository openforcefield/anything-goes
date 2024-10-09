
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


parser = argparse.ArgumentParser(description="Run Evaluator client")
parser.add_argument(
    "--port", "-p",
    type=int,
    default=8998,
)


def main(
    port: int = 8000,
    data_set_path: str = "small-test-set.json",
    force_field_path: str = "openff-2.2.0.offxml",
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
    )

import argparse
import logging
import pathlib
import sys
import time

import os

from openff.evaluator.backends.dask_kubernetes import DaskKubernetesExistingBackend
from openff.evaluator.backends.backends import ComputeResources, PodResources
from openff.toolkit.utils import OPENEYE_AVAILABLE
from openff.evaluator.server import EvaluatorServer
from openff.units import unit


logger = logging.getLogger(__name__)
logging.basicConfig(stream=sys.stdout, level=logging.DEBUG)

parser = argparse.ArgumentParser()
parser.add_argument("--cluster-name", type=str, default="evaluator-lw")
parser.add_argument("--namespace", type=str, default="openforcefield")
parser.add_argument("--storage-path", type=str, default="/evaluator-storage")
parser.add_argument("--port", type=int, default=8998)




if __name__ == "__main__":
    args = parser.parse_args()
    logger.info(f"OpenEye is available: {OPENEYE_AVAILABLE}")

    # change directory to storage path
    os.chdir(args.storage_path)

    working_directory = os.path.abspath(
        os.path.join(args.storage_path, "working-directory")
    )


    calculation_backend = DaskKubernetesExistingBackend(
        cluster_name=args.cluster_name,
        namespace=args.namespace,
        cluster_port=8786,
        resources_per_worker=PodResources(
            number_of_threads=1,
            memory_limit=8 * unit.gigabytes,
            ephemeral_storage_limit=20 * unit.gigabytes,
            number_of_gpus=1,
            preferred_gpu_toolkit=ComputeResources.GPUToolkit.CUDA,
        )
    )

    logger.info(f"Calculating with backend {calculation_backend}")
    with calculation_backend:
        evaluator_server = EvaluatorServer(
            calculation_backend,
            working_directory=working_directory,
            port=args.port,
            delete_working_files=False,
        )
        logger.info("Starting server")
        evaluator_server.start(asynchronous=False)
        

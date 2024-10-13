# Deploying Bespokefit as a webservice

This is a proof of concept method to deploy bespokefit as a webservice reachable via secure HTTPS connection managed by a traefik reverse proxy. 

We have used this at ASAPDiscovery to deploy inline forcefield generation for our Free Energy Calculation pipelines, using this setup on AWS


# IMPORTANT NOTE

Another reminder that this is provided under the terms of the MIT LICENSE as is, you are responsible for security on your own web deployments. Do not rely on me to have gotten everything right and audit carefully. We use HTTP to HTTPS redirection to try and ensure as much saftey as possible. 

# HOWTO

* copy the three files in this folder into the root of `openff-bespokefit
* you probably need to add the requisite drivers (eg torsiondirve) to the relevant environment file that docker builds from. I haven't had a chance to do this properly.
* edit the .env file and fill in the relevant variables
* run `docker-compose up -d` 
* check that the service is reachable on the client side with the below env file 


## Client side env file

```
# HTTPS default port
export BEFLOW_GATEWAY_PORT=443
export BEFLOW_GATEWAY_ADDRESS="https://asap-bespokefit.asapdata.org"
export BEFLOW_API_TOKEN="12345678910"
``

then run the below, checking that everything ran properly.

```
openff-bespoke executor  list
openff-bespoke executor submit --smiles CCCCCC --workflow debug
openff-bespoke executor watch --id 1
```



# AWS notes 

* Use an EC2 instance with at least 4 cores, 1 for the gateway and 3 for workers, 1 of each type. 
* Assign an elastic IP to the instance so it doesn't get re-jigged on restart

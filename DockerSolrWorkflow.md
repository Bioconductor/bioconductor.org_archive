Create a named container, e.g. "dev\_solr", as below : 

```
docker run -it makuk66/docker-solr:5.2.1 --name dev_solr
```
(Optionally) run some initial customization : 
```
docker exec -u root -it dev_solr /bin/bash -c "`cat DockerKickstart.sh`" /bin/bash -c "`cat DockerKickstart.sh`"
```
Login with the solr user account: 
```
docker exec -it dev_solr /bin/bash
```
Login as root:
```
docker exec -u root -it dev_solr /bin/bash
```
Login as root and get a real TTY (this is a workaround to a Docker issue) :
- https://github.com/docker/docker/issues/8755#issuecomment-83403289 

```
docker exec -u root -it gloomy_hodgkin script -q -c "/bin/bash" /dev/null
```

Create an example "core" (search index): 

```
# Command
docker exec -it --user=solr dev_solr bin/solr create_core -c gettingstarted

# Output
Setup new core instance directory:
/opt/solr/server/solr/gettingstarted

Creating new core 'gettingstarted' using command:
http://localhost:8983/solr/admin/cores?action=CREATE&name=gettingstarted&instanceDir=gettingstarted

{
  "responseHeader":{
    "status":0,
    "QTime":787},
  "core":"gettingstarted"}

```

Ingest some data : 
```
docker exec -it --user=solr my_solr bin/post -c gettingstarted example/exampledocs/manufacturers.xml
```

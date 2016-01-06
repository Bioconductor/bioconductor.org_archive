Create a named container, e.g. "dev\_solr", as below : 

```
docker run --name dev_solr -p 8983:8983 -it makuk66/docker-solr:5.2.1
```
(Optionally) run some initial customization : 
```
docker exec -u root -it dev_solr /bin/bash -c "`cat DockerKickstart.sh`"
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
docker exec -u root -it dev_solr script -q -c "/bin/bash" /dev/null
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
docker exec -it --user=solr dev_solr bin/post -c gettingstarted example/exampledocs/manufacturers.xml

docker exec -it --user=solr dev_solr bin/post -c gettingstarted example/exampledocs/sample.html
```

To ingest a large collection of data, do this : 
```
bin/post -c gettingstarted docs/
```

#### To query for data : 
- Open a browser, to the core's query page : 
    http://localhost:8983/solr/#/gettingstarted/query

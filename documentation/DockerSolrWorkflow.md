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


#### Testing ingest of production data, using "vanilla" defaults
1. Create a new core
```
docker exec -it --user=solr dev_solr bin/solr create_core -c test1
```
2. Ingest some files (this doesn't really work...)
```
curl -s "http://localhost:8983/solr/test1/update/extract" -F "myfile=@extra/www/bioc/packages/release/extra/aCGH.html"
curl -s "http://localhost:8983/solr/test1/update/extract" -F "myfile=@extra/www/bioc/packages/release/extra/ag.html"
curl -s "http://localhost:8983/solr/test1/update/extract" -F "myfile=@extra/www/bioc/packages/release/extra/adme16cod.html"
curl -s "http://localhost:8983/solr/test1/update/extract" -F "myfile=@extra/www/bioc/packages/release/extra/ath1121501.html"
```


1. Create a new core
```
docker exec -it --user=solr dev_solr bin/solr create_core -c test2
```
2. Ingest some files (results in decent output, but no path to file)
```
curl "http://localhost:8983/solr/test2/update/extract?commit=true&captureAttr=true" -F "myfile=@extra/www/bioc/packages/release/extra/html/ncdf.html"
curl "http://localhost:8983/solr/test2/update/extract?commit=true&captureAttr=true" -F "myfile=@extra/www/bioc/packages/release/extra/html/XMLRPC.html"
curl "http://localhost:8983/solr/test2/update/extract?commit=true&captureAttr=true" -F "myfile=@extra/www/bioc/packages/release/extra/html/RCurl.html""

```


1. Create a new core
```
docker exec -it --user=solr dev_solr bin/solr create_core -c test3
```
2. Ingest some files 
```
curl "http://localhost:8983/solr/test3/update/extract?commit=true&captureAttr=true&resource.name=extra/www/bioc/packages/release/extra/html/ncdf.html" -F "myfile=@extra/www/bioc/packages/release/extra/html/ncdf.html"
curl "http://localhost:8983/solr/test3/update/extract?commit=true&captureAttr=true&resource.name=extra/www/bioc/packages/release/extra/html/XMLRPC.html" -F "myfile=@extra/www/bioc/packages/release/extra/html/XMLRPC.html"
curl "http://localhost:8983/solr/test3/update/extract?commit=true&captureAttr=true&resource.name=extra/www/bioc/packages/release/extra/html/RCurl.html" -F "myfile=@extra/www/bioc/packages/release/extra/html/RCurl.html"
```
3. The above is sufficient, to allow us to execute a query like : 
`http://localhost:8983/solr/test3/select?q=*%3A*&fl=resourcename%2Ctitle%2Ccontent&wt=json&indent=true`
and retrieve results like :
```
{
  "responseHeader":{
    "status":0,
    "QTime":0,
    "params":{
      "q":"*:*",
      "indent":"true",
      "fl":"resourcename,title,content",
      "wt":"json"}},
  "response":{"numFound":3,"start":0,"docs":[
      {
        "resourcename":["extra/www/bioc/packages/release/extra/html/ncdf.html"],
        "title":["ncdf"]},
      {
        "resourcename":["extra/www/bioc/packages/release/extra/html/XMLRPC.html"],
        "title":["XMLRPC"]},
      {
        "resourcename":["extra/www/bioc/packages/release/extra/html/RCurl.html"],
        "title":["RCurl"]}]
  }}
```




1. Create a new core
```
docker exec -it --user=solr dev_solr bin/solr create_core -c test4
```
2. Ingest some files 
```
curl "http://localhost:8983/solr/test4/update/extract?extractOnly=true&captureAttr=true&resource.name=extra/www/bioc/packages/release/extra/html/ncdf.html" -F "myfile=@extra/www/bioc/packages/release/extra/html/ncdf.html"
curl "http://localhost:8983/solr/test4/update/extract?extractOnly=true&captureAttr=true&resource.name=extra/www/bioc/packages/release/extra/html/XMLRPC.html" -F "myfile=@extra/www/bioc/packages/release/extra/html/XMLRPC.html"
curl "http://localhost:8983/solr/test4/update/extract?extractOnly=true&captureAttr=true&resource.name=extra/www/bioc/packages/release/extra/html/RCurl.html" -F "myfile=@extra/www/bioc/packages/release/extra/html/RCurl.html"
```
3. The above is sufficient, to allow us to execute a query like : 
`http://localhost:8983/solr/test3/select?q=*%3A*&fl=resourcename%2Ctitle%2Ccontent&wt=json&indent=true`
and retrieve results like :
```
{
  "responseHeader":{
    "status":0,
    "QTime":0,
    "params":{
      "q":"*:*",
      "indent":"true",
      "fl":"resourcename,title,content",
      "wt":"json"}},
  "response":{"numFound":3,"start":0,"docs":[
      {
        "resourcename":["extra/www/bioc/packages/release/extra/html/ncdf.html"],
        "title":["ncdf"]},
      {
        "resourcename":["extra/www/bioc/packages/release/extra/html/XMLRPC.html"],
        "title":["XMLRPC"]},
      {
        "resourcename":["extra/www/bioc/packages/release/extra/html/RCurl.html"],
        "title":["RCurl"]}]
  }}
```

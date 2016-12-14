# Solr Backup

For the purposes of recovery & development, we need a backup of our Solr data. [1.]
One way to backup the Solr data is described below.

Backup:
```
# SSH to the master node, and forward port 8983 from the local machine to master
ssh webadmin@master.bioconductor.org -L8983:localhost:8983

# Create a backup directory
mkdir ~/SOLR-BACKUP-DIR

# Invoke the backup endpoint.  This will create a backup
# of the 'default' core at '/home/webadmin/SOLR-BACKUP-DIR'.
curl 'http://localhost:8983/solr/default/replication?command=backup&location=/home/webadmin/SOLR-BACKUP-DIR'

# View the status of the backup
curl 'http://localhost:8983/solr/default/replication?command=details&indent=on'
```
_In my test, this creates a directory `/home/webadmin/SOLR-BACKUP-DIR/snapshot.20160113195440762` on master._


Restoration:
```
# We'll restore the 'default' core to a new core
# called 'restore-test'.  To do so, we need to create it:
$HOME/solr-5.2.1/bin/solr create -c restore-test
```
If that succeeds, you should see output like this :
```
Setup new core instance directory:
/home/webadmin/solr-5.2.1/server/solr/restore-test

Creating new core 'restore-test' using command:
http://localhost:8983/solr/admin/cores?action=CREATE&name=restore-test&instanceDir=restore-test

{
  "responseHeader":{
    "status":0,
    "QTime":470},
  "core":"restore-test"}

```
Next, you'll restore the previous backup (`snapshot.20160113195440762`) to our new
core (`restore-test`) :
```
curl 'http://localhost:8983/solr/restore-test/replication?command=restore&name=snapshot.20160113195440762&indent=on'
```

# FIXME
Browsing to the web interface, there are no documents : http://localhost:8983/solr/#/restore-test .

Invoking commit made no difference :
```
curl http://localhost:8983/solr/restore-test/update --data-binary '<commit/>' -H 'Content-type:text/xml; charset=utf-8'
```
Refresh the Solr web interface and see that there are still no documents.

[^1.]: Based on Solr backup documentation: https://cwiki.apache.org/confluence/display/solr/Making+and+Restoring+Backups+of+SolrCores

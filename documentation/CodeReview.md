
#### Code Review 2016-Jan-07
The following notes were captured by a code review performed by Martin and Brian.  We
need to confirm the below is correct with Dan.  This information relates to the
[Site Search](../README.md#site-search).  Perhaps once we discuss this and confirm our
understanding, we can integrate the following into the current "_Site Search_"
section of the [README](../README.md#site-search).

##### Our understanding of the code
- Indexing is kicked off by a cron job that's on staging.
- Cron job runs 'index_production' task from the Rakefile.
- 'index_production' syncs scripts to master
- 'index_production' runs get_links.rb and 'do_index.rb'
- 'do_index.rb' isn't in source control, and it runs get_links.rb itself.
- 'do_index.rb' sources search_indexer.rb, and constructs it (the constructor does work)
- search_indexer's constructor: creates "index.sh" with curl commands to interact with
  Solr (over HTTP)
- A file named search_indexer_cache.yaml exists, and it seems to have k/v pairs
  identifying things that are known to Solr.
  - This is called "cache" in search_indexer.rb and we think it might be used to improve
    performance (short circuit) the process of performing the day's search index update.  
    The update, which is done by index.sh, currently has ~40,000 invocations of curl,
    which might indicate the caching mechanism isn't working properly.
- 'search_indexer.rb' has a function 'get_boost' (meant to improve score of release +100
   and devel +50)
- cachecopy in 'search_indexer.rb' is a hash of booleans, indicating which items in the
  cache (another hash) are meant to be indexed this is based upon their file extension.
- 'search_indexer.rb' determines which files are to be deleted from Solr, by subsetting :
  `to_be_deleted = cache.keys.sort - cachecopy.keys.sort`.
- 'search_indexer.rb' finishes index.sh with a <commit/> line (curl, communicating with
  Solr), and makes index.sh executable.
- Next, 'index_production' (in the Rakefile) invokes index.sh, which was created.
- 'index_production' (in the Rakefile) is done.

##### Improvements to be made at some point
- Cleanup: Remove old code that is commented out
- Cleanup: Remove old functions that aren't used (e.g. `get_list_of_files_to_index_old`
  in 'search_indexer.rb'
- Cleanup: Why are we creating a shell script in 'search_indexer.rb', then invoking it
    from the Rakefile?  This is a single unit, should be a single function.

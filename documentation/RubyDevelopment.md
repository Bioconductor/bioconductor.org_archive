Ruby Development Steps
======================
Ensure rbenv is installed : 
```
# *nix
which rbenv
# Windows
where rbenv
```
Install Ruby:
- I'm using version 2.1.1 using rbenv.
- Production seems to use ruby 1.9.3p484 (2013-11-22 revision 43786) [x86_64-linux]

NB: The README incorrectly says that Ruby the site requires ruby 1.9.1 or newer.  However,
the dependency on `json -v '1.8.1'` will cause `bundle install` to fail when building on
Ruby 2.3.0 .

```
blong@work:~/Documents/Work/REPOS__git/SOLR/bioconductor.org$ rbenv install 2.3.0
Downloading ruby-2.3.0.tar.bz2...
-> https://cache.ruby-lang.org/pub/ruby/2.3/ruby-2.3.0.tar.bz2
Installing ruby-2.3.0...
Installed ruby-2.3.0 to /home/blong/.rbenv/versions/2.3.0
```

Set the current version as global 
```
blong@work:~/Documents/Work/REPOS__git/SOLR/bioconductor.org$ rbenv global 2.3.0
```
Install bundler
```
blong@work:~/Documents/Work/REPOS__git/SOLR/bioconductor.org$ gem install bundler
Fetching: bundler-1.11.2.gem (100%)
Successfully installed bundler-1.11.2
Parsing documentation for bundler-1.11.2
Installing ri documentation for bundler-1.11.2
Done installing documentation for bundler after 2 seconds
1 gem installed
```
Install bioconductor.org 
```
blong@work:~/Documents/Work/REPOS__git/SOLR/bioconductor.org$ bundle install
Fetching gem metadata from http://rubygems.org/..........
Fetching version metadata from http://rubygems.org/...
Fetching dependency metadata from http://rubygems.org/..
Using rake 10.4.2
Installing addressable 2.3.6
Installing rack 1.6.0
Installing algorithms 0.5.0 with native extensions
Installing json 1.8.1 with native extensions

Gem::Ext::BuildError: ERROR: Failed to build gem native extension.

    current directory: /home/blong/.rbenv/versions/2.3.0/lib/ruby/gems/2.3.0/gems/json-1.8.1/ext/json/ext/generator
/home/blong/.rbenv/versions/2.3.0/bin/ruby -r ./siteconf20160104-24891-q8lz1w.rb extconf.rb
creating Makefile

current directory: /home/blong/.rbenv/versions/2.3.0/lib/ruby/gems/2.3.0/gems/json-1.8.1/ext/json/ext/generator
make "DESTDIR=" clean

current directory: /home/blong/.rbenv/versions/2.3.0/lib/ruby/gems/2.3.0/gems/json-1.8.1/ext/json/ext/generator
make "DESTDIR="
compiling generator.c
In file included from generator.c:1:0:
../fbuffer/fbuffer.h: In function ‘fbuffer_to_s’:
../fbuffer/fbuffer.h:175:47: error: macro "rb_str_new" requires 2 arguments, but only 1 given
     VALUE result = rb_str_new(FBUFFER_PAIR(fb));
                                               ^
../fbuffer/fbuffer.h:175:20: warning: initialization makes integer from pointer without a cast [enabled by default]
     VALUE result = rb_str_new(FBUFFER_PAIR(fb));
                    ^
make: *** [generator.o] Error 1

make failed, exit code 2

Gem files will remain installed in /home/blong/.rbenv/versions/2.3.0/lib/ruby/gems/2.3.0/gems/json-1.8.1 for inspection.
Results logged to /home/blong/.rbenv/versions/2.3.0/lib/ruby/gems/2.3.0/extensions/x86_64-linux/2.3.0-static/json-1.8.1/gem_make.out
Installing mini_portile 0.6.2
Installing buftok 0.2.0
Installing columnize 0.9.0
Installing coderay 1.1.0
Installing colored 1.2
Installing descriptive_statistics 2.5.1
Installing unf_ext 0.0.6 with native extensions
Installing equalizer 0.0.9
Installing multipart-post 2.0.0
Installing tilt 2.0.1
Installing hpricot 0.8.6 with native extensions
Installing htmlentities 4.3.3
Installing http_parser.rb 0.6.0 with native extensions
Installing multi_xml 0.5.5
Installing kramdown 1.6.0
Installing systemu 2.6.4
Installing mime-types 2.4.3
Installing net-http-digest_auth 1.4
Installing net-http-persistent 2.9.4
Installing ntlm-http 0.1.1
Installing webrobots 0.1.1
Installing thread_safe 0.3.4
Installing method_source 0.8.2
Installing naught 1.0.0
Installing pg 0.18.1 with native extensions
Installing polyglot 0.3.5
Installing slop 3.6.0
Installing rdiscount 2.1.7.1 with native extensions
Installing redis 3.2.0
Installing stream 0.5
Installing sass 3.4.9
Installing sequel 4.22.0
Installing simple_oauth 0.3.1
Installing sqlite3 1.3.10 with native extensions
Using bundler 1.11.2
Installing adsf 1.2.0
Installing rack-cache 1.2
An error occurred while installing json (1.8.1), and Bundler cannot continue.
Make sure that `gem install json -v '1.8.1'` succeeds before bundling.

```

#!/usr/bin/env bash

if [ "$EUID" -ne 0 ]; then
  echo "Please run as root"
  exit 1
else
  cat <<'EOF'
    
  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  Good news
  
  We're running DockerKickstart as root!
  
  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
EOF

fi

apt-get -o Acquire::ForceIPv4=true update
apt-get -o Acquire::ForceIPv4=true install -y vim

### Customize Bash
# Enable colorized output, by adding this to
# the bottom of ... uh... somewhere.
# FIXME: Find the right file!   I already tried 
#    - /etc/bash.bashrc
#    - /etc/profile
#    - /etc/profile.d/docker-kickstart.sh

LS_OPTIONS='--color=auto'

cat << EOF >> /etc/bash.bashrc

alias ls='ls $LS_OPTIONS'
alias ll='ls $LS_OPTIONS -l'
alias l='ls $LS_OPTIONS -lA'
EOF


### Customize Vim
# Enable syntax highlighting
echo 'syntax on' >> /etc/vim/vimrc

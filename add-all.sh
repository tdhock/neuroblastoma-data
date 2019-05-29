#!/bin/bash
git status|grep data|grep -v new|grep -v modified|xargs git add

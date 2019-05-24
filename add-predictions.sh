#!/bin/bash
git status|grep data|grep -v new|grep -v modified|sed 's#models/#models/*/predictions.csv#'|xargs git add

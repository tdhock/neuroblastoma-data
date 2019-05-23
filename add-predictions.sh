#!/bin/bash
git status|grep data|grep -v new|sed 's#models/#models/*/predictions.csv#'|xargs git add

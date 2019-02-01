#!/usr/bin/env bash

prodigy ner.manual mindfields_ner en_core_web_sm raw_corpus.jsonl --label "ORGANISM, ACTIVITY, CONCENTRATION, CONCENTRATION_TYPE, ACTIVITY_MODULATION"
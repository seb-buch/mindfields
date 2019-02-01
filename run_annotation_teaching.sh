#!/usr/bin/env bash
prodigy ner.teach mindfields_ner model raw_corpus.jsonl --label "ORGANISM, ACTIVITY, CONCENTRATION, CONCENTRATION_TYPE, ACTIVITY_MODULATION"
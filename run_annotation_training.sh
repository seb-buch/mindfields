#!/usr/bin/env bash
prodigy ner.batch-train mindfields_ner en_core_web_sm --output ./model --eval-split 0.5 --label "ORGANISM, ACTIVITY, CONCENTRATION, CONCENTRATION_TYPE, ACTIVITY_MODULATION"
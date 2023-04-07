#!/bin/bash

dir="/Users/badgerjh/lib/src/llama.cpp"
main="$dir/main"
model="$dir/models/ggml-vicuna-13b-4bit-rev1.bin"
prompt="$dir/prompts/chat-with-vicuna.txt"
revprompt="### Human:"
temp=0.8
n=-1
k=40

POSITIONAL_ARGS=()

while [[ $# -gt 0 ]]; do
  case $1 in
    -m|--model)
      model="$2"
      shift # past argument
      shift # past value
      ;;
    -n|--number_of_tokens)
      n=$2
      shift # past argument
      shift # past value
      ;;
    -k|--top_k)
      k=$2
      shift # past argument
      shift # past value
      ;;
    -p|--prompt)
      prompt="$2"
      shift # past argument
      shift # past value
      ;;
    -r|--reverse_prompt)
      revprompt="$2"
      shift # past argument
      shift # past value
      ;;
    -t|--temp)
      temp=$2
      shift # past argument
      shift # past value
      ;;
    *)
      POSITIONAL_ARGS+=("$1") # save positional arg
      shift # past argument
      ;;
  esac
done



cmd="rlwrap $main -m $model --color -f $prompt -i -r '$revprompt' --top_k $k -t 8 --temp $temp --ignore-eos -n $n -c 2048 -t 7"
eval $cmd

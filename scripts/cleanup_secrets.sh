#!/usr/bin/sh
# shellcheck disable=SC3040
set -euo pipefail

usage() {
cat << EOF
Usage: $0 [-a <aws profile> -r <region>]

Deletes all secrets with 7 day recovery window for secrets that match "tower-.*\/TBCHASIN_.*".
NB: Seqera Platform 24.2.0 prefixes secrets with "tower-<job_id>/" in AWS secretsmanager.

EOF
}

region='us-east-1'

while getopts "a:r:" OPTION; do
    case $OPTION in
        a) aws_profile=$OPTARG;;
        r) region=$OPTARG;;
        ?) usage; exit 1;
    esac
done

secret_names=$(aws ${aws_profile:+--profile $aws_profile} ${region:+--region $region} secretsmanager list-secrets | jq -r '.[] | .[] | .Name')

for secret_name in $secret_names; do
    # shellcheck disable=SC3010
    if [[ "$secret_name" =~ tower-.*\/TBCHASIN_.* ]]; then

        # shellcheck disable=SC3045
        # shellcheck disable=SC2162
        read -n 1 -p "Confirm delete $secret_name [y/n]: " confirm

        echo

        if [ "$confirm" = "y" ]; then
            echo "Deleting $secret_name"
            aws ${aws_profile:+--profile $aws_profile} ${region:+--region $region} \
                secretsmanager delete-secret \
                --secret-id "$secret_name" \
                --recovery-window-in-days 7
        else
            echo "Did not delete $secret_name"
        fi
    fi
done
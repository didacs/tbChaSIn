# AWS setup for ECR respositories

Deploy and initialize ECR images with:

```console
aws [--profile <profile>] cloudformation deploy --stack-name tbChaSIn-ECR-repositories --template aws/tbChaSIn-resources.yml --region us-east-1
```

Update the Cloudformation stack by using the AWS console and uploading the updated CFN template.

> NB: Unfortunately, the `aws cloudformation update-stack` command can only take template URLs in S3.